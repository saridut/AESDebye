#include <IntensityKernelsGPU.hpp>

__global__ void calculateIntensityKernelGPU(int nBins, int nqSteps,
                 double *d_centers, double *d_counts, double *d_qVector, double *d_intensity)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nqSteps)
    {
        double q = d_qVector[i];

        double irq = 1.0;
        double rq = 0.0;
        double intensity = 0.0;
        for (uint binId=0; binId < nBins; binId++)
        {
            auto count = d_counts[binId];
            if (count == 0)
                continue;
            rq = d_centers[binId] * q;

            if (rq > FLT_EPSILON)
            {
                irq = sin(rq) / rq;
            }

            intensity += count * irq;
        }
        d_intensity[i] = intensity;
    }
}

std::vector<double>
calculateIntensityGPU(std::vector<double> const& qVector,
                      std::vector<double> const& centers,
                      const std::vector<double> &counts)
{

    // managed centers and counts
    double *d_centers;
    double *d_counts;
    double *d_intensity;
    double *d_qVector;
    size_t nqSteps = qVector.size();
    cudaMalloc(&d_centers, centers.size() * sizeof(double));
    cudaMalloc(&d_counts, counts.size() * sizeof(double));
    cudaMalloc(&d_intensity, nqSteps * sizeof(double));
    cudaMalloc(&d_qVector, nqSteps * sizeof(double));

    // copy centers and counts to device
    cudaMemcpy(d_centers, centers.data(), centers.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_counts, counts.data(), counts.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_qVector, qVector.data(), nqSteps * sizeof(double), cudaMemcpyHostToDevice);

    // calculate intensity
    dim3 threadsPerBlock(512);
    dim3 numBlocks((nqSteps + threadsPerBlock.x - 1) / threadsPerBlock.x);

    calculateIntensityKernelGPU<<<numBlocks, threadsPerBlock>>>(centers.size(), nqSteps,
            d_centers,d_counts, d_qVector, d_intensity);

    cudaDeviceSynchronize();

    // copy intensity to host
    std::vector<double> intensity(nqSteps);
    cudaMemcpy(intensity.data(), d_intensity, nqSteps * sizeof(double), cudaMemcpyDeviceToHost);

    // free memory
    cudaFree(d_centers);
    cudaFree(d_counts);
    cudaFree(d_intensity);
    cudaFree(d_qVector);

    return intensity;
}