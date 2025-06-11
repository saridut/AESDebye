#include <IntensityKernels.hpp>

#pragma region reduction

#pragma omp declare                                                                                                                         \
    reduction(                                                                                                                              \
            + : std::vector<double> : std::transform(omp_out.begin(),\
            omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

#pragma endregion reduction

#pragma region intensity_calculation

std::vector<double>
calculateIntensityCPU(std::vector<double> const& qVector,
                      std::vector<double> const& centers,
                      const std::vector<double> &counts)
{
    std::vector<double> intensity(qVector.size(), 0.0);
    uint nqSteps = qVector.size();
#pragma omp parallel for schedule(dynamic) default(none) shared(centers, counts, qVector, nqSteps) reduction(+ : intensity)
    for (uint binId = 0; binId < counts.size(); binId++)
    {
        if (counts[binId] == 0)
            continue;

        double count = counts[binId];
        double correctedCenter = centers[binId];

        double q, rq, irq;
        for (uint qIdx = 0; qIdx < nqSteps; qIdx++)
        {
            q = qVector[qIdx];
            irq = 1.0;
            if (q > FLT_EPSILON)
            {
                rq = correctedCenter * q;
                if (rq > FLT_EPSILON)
                    irq = sin(rq) / rq;
            }
            intensity[qIdx] += (irq * count);
        }
    }
    return intensity;
}

#pragma endregion intensity_calculation