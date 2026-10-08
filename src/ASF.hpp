#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>

/**
 * @brief Supported formulation standards for Atomic Scattering Factor (ASF) calculations.
 * - WaasmaierKirfel5: 5 Gaussians + constant (Waasmaier & Kirfel 1995, standard in DebUsSy).
 *                     Formula: f0(q) = c + SUM_{i=1..5} a_i * exp(-b_i * (q / 4pi)^2)
 * - CromerMann4:      4 Gaussians + constant (Cromer & Mann 1968 / International Tables Vol C).
 *                     Formula: f0(q) = c + SUM_{i=1..4} a_i * exp(-b_i * (q / 4pi)^2)
 */
enum class ASFFormulation {
    WaasmaierKirfel5,
    CromerMann4
};

inline std::string asfFormulationToString(ASFFormulation form) {
    switch (form) {
        case ASFFormulation::WaasmaierKirfel5:
            return "WaasmaierKirfel5";
        case ASFFormulation::CromerMann4:
            return "CromerMann4";
    }
    return "WaasmaierKirfel5";
}

inline ASFFormulation stringToASFFormulation(const std::string& str) {
    std::string lower = str;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "cromermann4" || lower == "cm4" || lower == "cm" || lower == "4") {
        return ASFFormulation::CromerMann4;
    }
    // Default is WaasmaierKirfel5
    return ASFFormulation::WaasmaierKirfel5;
}

/**
 * @brief Struct holding analytical coefficients for an element's atomic form factor.
 * Supports up to 5 Gaussian components (a1..a5, b1..b5) and constant c.
 * Also holds optional metadata: neutron bound coherent scattering length, Z, and atomic radius.
 */
struct ASFCoeffs {
    double a1{0}, b1{0}, a2{0}, b2{0}, a3{0}, b3{0}, a4{0}, b4{0}, a5{0}, b5{0}, c{0};
    double coh_b{0};          // Bound coherent neutron scattering length [fm]
    int atomic_number{0};     // Z
    double atomic_radius{0};  // Radius in Angstrom

    ASFCoeffs() = default;

    // 4-Gaussian constructor (Cromer-Mann)
    ASFCoeffs(double a1, double b1, double a2, double b2, double a3, double b3, double a4, double b4, double c)
        : a1(a1), b1(b1), a2(a2), b2(b2), a3(a3), b3(b3), a4(a4), b4(b4), a5(0), b5(0), c(c) {}

    // 5-Gaussian constructor (Waasmaier-Kirfel)
    ASFCoeffs(double a1, double a2, double a3, double a4, double a5, double c,
              double b1, double b2, double b3, double b4, double b5,
              double coh_b = 0.0, int atomic_number = 0, double atomic_radius = 0.0)
        : a1(a1), b1(b1), a2(a2), b2(b2), a3(a3), b3(b3), a4(a4), b4(b4), a5(a5), b5(b5), c(c),
          coh_b(coh_b), atomic_number(atomic_number), atomic_radius(atomic_radius) {}

    [[nodiscard]] inline double ASFq(double q) const
    {
        // Factor = (sin(theta)/lambda)^2 = (q / (4*pi))^2 = q^2 / (16 * pi^2)
        const double factor = q * q * 0.06250 * M_1_PI * M_1_PI;
        double asf = c;
        asf += a1 * std::exp(-b1 * factor);
        asf += a2 * std::exp(-b2 * factor);
        asf += a3 * std::exp(-b3 * factor);
        asf += a4 * std::exp(-b4 * factor);
        if (a5 != 0.0 || b5 != 0.0)
        {
            asf += a5 * std::exp(-b5 * factor);
        }
        return asf;
    }
};

/**
 * @brief Table containing atomic scattering factor coefficients for elements.
 * Loaded cleanly from a TSV/CSV/whitespace-delimited file.
 */
class ASFTable {
private:
    std::unordered_map<std::string, ASFCoeffs> table;
    ASFFormulation formulation{ASFFormulation::WaasmaierKirfel5};
    std::string loadedFilePath{""};

public:
    static constexpr const char* DEFAULT_ASF_FILE = "src/aesdebye/data/asf.tsv";

    ASFTable() = default;

    explicit ASFTable(const std::string& filepath, ASFFormulation form = ASFFormulation::WaasmaierKirfel5)
    {
        loadFromFile(filepath.empty() ? DEFAULT_ASF_FILE : filepath, form);
    }

    /**
     * @brief Load ASF table from a TSV, CSV, or whitespace-delimited file.
     * Enforces the specified formulation (WaasmaierKirfel5 or CromerMann4).
     *
     * @param filepath Path to the ASF file.
     * @param form Formulation to parse and enforce (default: WaasmaierKirfel5).
     * @throws std::runtime_error on file open error or parsing failure.
     */
    void loadFromFile(const std::string& filepath, ASFFormulation form = ASFFormulation::WaasmaierKirfel5)
    {
        std::ifstream file(filepath);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open ASF file: " + filepath);
        }

        std::unordered_map<std::string, ASFCoeffs> newTable;
        std::string line;
        size_t lineNum = 0;

        try
        {
            while (std::getline(file, line))
            {
                lineNum++;
                // Skip empty lines and comments
                if (line.empty()) continue;
                size_t firstNonBlank = line.find_first_not_of(" \t\r\n");
                if (firstNonBlank == std::string::npos || line[firstNonBlank] == '#')
                {
                    continue;
                }

                // Clean delimiters (commas, brackets, colons)
                std::string processedLine = line;
                for (char &ch : processedLine)
                {
                    if (ch == ',' || ch == '[' || ch == ']' || ch == ':')
                    {
                        ch = ' ';
                    }
                }

                std::istringstream iss(processedLine);
                std::string element;
                if (!(iss >> element)) continue;

                std::vector<std::string> tokens;
                std::string token;
                while (iss >> token)
                {
                    tokens.push_back(token);
                }

                auto parseVal = [](const std::string& s) -> double {
                    if (s == "nan" || s == "NAN" || s == "NaN" || s == "null" || s == "NULL") {
                        return 0.0;
                    }
                    return std::stod(s);
                };

                if (form == ASFFormulation::WaasmaierKirfel5)
                {
                    // Expect at least 11 values: a1..a5, c, b1..b5
                    if (tokens.size() < 11)
                    {
                        throw std::runtime_error("Line " + std::to_string(lineNum) +
                                                 ": Expected at least 11 values for WaasmaierKirfel5 (a1..a5 c b1..b5), found " +
                                                 std::to_string(tokens.size()) + " for element " + element);
                    }

                    double a1 = parseVal(tokens[0]);
                    double a2 = parseVal(tokens[1]);
                    double a3 = parseVal(tokens[2]);
                    double a4 = parseVal(tokens[3]);
                    double a5 = parseVal(tokens[4]);
                    double c  = parseVal(tokens[5]);
                    double b1 = parseVal(tokens[6]);
                    double b2 = parseVal(tokens[7]);
                    double b3 = parseVal(tokens[8]);
                    double b4 = parseVal(tokens[9]);
                    double b5 = parseVal(tokens[10]);

                    double coh_b = 0.0;
                    int atomic_number = 0;
                    double atomic_radius = 0.0;

                    if (tokens.size() > 11) coh_b = parseVal(tokens[11]);
                    if (tokens.size() > 12) atomic_number = static_cast<int>(parseVal(tokens[12]));
                    if (tokens.size() > 13) atomic_radius = parseVal(tokens[13]);

                    newTable[element] = ASFCoeffs(a1, a2, a3, a4, a5, c, b1, b2, b3, b4, b5, coh_b, atomic_number, atomic_radius);
                }
                else // CromerMann4
                {
                    // Expect 9 values: a1, b1, a2, b2, a3, b3, a4, b4, c
                    if (tokens.size() < 9)
                    {
                        throw std::runtime_error("Line " + std::to_string(lineNum) +
                                                 ": Expected at least 9 values for CromerMann4 (a1 b1 a2 b2 a3 b3 a4 b4 c), found " +
                                                 std::to_string(tokens.size()) + " for element " + element);
                    }

                    double a1 = parseVal(tokens[0]);
                    double b1 = parseVal(tokens[1]);
                    double a2 = parseVal(tokens[2]);
                    double b2 = parseVal(tokens[3]);
                    double a3 = parseVal(tokens[4]);
                    double b3 = parseVal(tokens[5]);
                    double a4 = parseVal(tokens[6]);
                    double b4 = parseVal(tokens[7]);
                    double c  = parseVal(tokens[8]);

                    newTable[element] = ASFCoeffs(a1, b1, a2, b2, a3, b3, a4, b4, c);
                }
            }
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error("Error parsing ASF file '" + filepath + "': " + e.what());
        }

        if (newTable.empty())
        {
            throw std::runtime_error("No valid element entries found in ASF file: " + filepath);
        }

        table = std::move(newTable);
        formulation = form;
        loadedFilePath = filepath;
    }

    [[nodiscard]] bool contains(const std::string& element) const
    {
        if (element.empty() || element == "None" || element == "none" || element == "NONE")
        {
            return true;
        }
        return table.find(element) != table.end();
    }

    [[nodiscard]] const ASFCoeffs& get(const std::string& element) const
    {
        if (element.empty() || element == "None" || element == "none" || element == "NONE")
        {
            static const ASFCoeffs unityCoeffs{};
            return unityCoeffs;
        }
        auto it = table.find(element);
        if (it != table.end())
        {
            return it->second;
        }
        throw std::invalid_argument("Element '" + element + "' not found in ASF table (loaded from: " +
                                    (loadedFilePath.empty() ? DEFAULT_ASF_FILE : loadedFilePath) + ")");
    }

    [[nodiscard]] const ASFCoeffs& operator[](const std::string& element) const
    {
        return get(element);
    }

    [[nodiscard]] size_t size() const { return table.size(); }
    [[nodiscard]] ASFFormulation getFormulation() const { return formulation; }
    [[nodiscard]] std::string getLoadedFilePath() const { return loadedFilePath; }
};
