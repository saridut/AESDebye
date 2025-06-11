#include <Profile.hpp>
#include <asciichart/include/ascii/ascii.h>


Profile Profile::operator+(const Profile &rhs) const
{
    // create new intensity
    std::vector<double> newI = intensity;

    // add rhs intensity to new intensity
    #pragma omp parallel for schedule(static) default(none) shared(newI, rhs)
    for (size_t i = 0; i < newI.size(); i++)
    {
        newI[i] += rhs.intensity[i];
    }

    // return new profile
    return {q, newI, twoTheta, calculationTime + rhs.calculationTime, testPassed && rhs.testPassed};
}

Profile Profile::operator-(const Profile &rhs) const
{       
        // new intensity
        std::vector<double> newI = intensity;

        // subtract rhs intensity from new intensity
        #pragma omp parallel for schedule(static) default(none) shared(newI, rhs)
        for (size_t i = 0; i < newI.size(); i++)
        {
            newI[i] -= rhs.intensity[i];
        }

        // return new profile
        return {q, newI,twoTheta,calculationTime + rhs.calculationTime, testPassed && rhs.testPassed};
}


bool Profile::test()
{   
    int64 sum = 0;
    for (auto i : intensity)
    {
        if (i < 0)     // check if intensity is negative
        {
            testPassed = false;
            return testPassed;
        }
        sum += i;
    }
    testPassed = sum > 0; // check if the total intensity > 0;
    return testPassed;
}


void Profile::toCSV(std::string filename)
{
    std::ofstream file;
    file.open(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    // header 
    file << "q,2theta,intensity" << std::endl;
    for (int i = 0; i < q.size(); i++)
    {
        file << q[i] << "," << twoTheta[i] << "," << intensity[i] << std::endl;
    }
    file.close();
}


std::string Profile::termPlot(int skip, bool log) const
{
    // empty container
    std::vector<double> series;
    double total = 0;
    skip = intensity.size() / 100;
    for (int i = 0; i < intensity.size(); i += skip)
    {
        // add the log of the intensity
        double I = intensity.at(i) > 0 ? intensity.at(i) : 1e-0;
        series.push_back(log ? std::log10(I) : I);
    }

    // create the ascii chart
    ascii::Asciichart asciichart(std::vector<std::vector<double>>{series});

    // plot it on the first rank
    return asciichart.type(ascii::Asciichart::LINE).height(10).Plot();
}