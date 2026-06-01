#include <utility>

#include <DataTypes.hpp>


class Profile
{
public:
    std::vector<double> q={};
    std::vector<double> intensity={};
    std::vector<double> twoTheta={};
    double calculationTime=0.0;
    bool testPassed=false;

    Profile() = default;
    Profile(std::vector<double> q,
            std::vector<double> intensity,
            std::vector<double> twoTheta,
            double calculationTime,
            bool testPassed) : q(std::move(q)), intensity(std::move(intensity)), twoTheta(std::move(twoTheta)), calculationTime(calculationTime), testPassed(testPassed) {}

    Profile(double start, double end, size_t nSteps, bool twoThetaSpace, double lambda)
    {
        q = std::vector<double>(nSteps);
        intensity = std::vector<double>(nSteps);
        twoTheta = std::vector<double>(nSteps);
        double stepSize = (end - start) / (double) nSteps;
        bool warn = false;

#pragma omp parallel for schedule(static) shared(q, twoTheta, twoThetaSpace, start, end, nSteps, stepSize, lambda, warn) default(none)
        for (int i = 0; i < nSteps; i++)
        {
            if (twoThetaSpace)
            {
                double radians = (start + i * stepSize) * M_PI / 180;
                q[i] = 4 * M_PI * std::sin(radians / 2) / lambda;
                twoTheta[i] = start + i * stepSize;
            }
            else
            {
                q[i] = start + i * stepSize;
                double value = q[i] * lambda / 4 / M_PI;
                if (value > 1.0 or value < -1.0)
                {
                    warn = true;
                    value = std::max(-1.0, std::min(1.0, value));
                }
                value = std::asin(value);
                twoTheta[i] = 2 * value * 180 / M_PI;
            }
        }

        if (warn)
        {
            std::cerr << "Warning: During the conversion from q to theta, some values are outside the valid range for the given lambda. "
                      << "Some values have been clipped to the valid range." << std::endl;
        }
    }

    Profile operator+(const Profile &rhs) const;
    Profile operator-(const Profile &rhs) const;

    std::string toString() const
    {
        std::stringstream ss;
        ss << "Test: " << helpers::bool_to_string(testPassed) << " time taken: " << calculationTime << "s" << "\n";
        ss << termPlot();
        return ss.str();
    }

    bool test();
    void toCSV(std::string filename);    
    [[nodiscard]] std::string termPlot(int skip=10, bool log=true) const;
};
