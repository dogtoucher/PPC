#include<cmath>
#include<iostream>
void correlate(int ny, int nx, const float* data, float* result) {
    double* mean = new double[ny]{0.0};
    double* std = new double[ny]{0.0};
    double* normalizedData = new double[ny*nx]{0.0};
    for (int i = 0; i < ny; i++) {
        for (int k = 0; k < nx; k++) {
            mean[i] += data[i*nx + k];
        }
        mean[i] /= nx;
    }
    for (int i = 0; i < ny; i++) {
        double sumSquares = 0.0;
        for (int k = 0; k < nx; k++) {
            sumSquares += pow(data[i*nx + k]-mean[i],2);
        }
        std[i] = sqrt(sumSquares);
    }
    for (int i = 0; i < ny; i++) {
        for (int k = 0; k < nx; k++) {
            normalizedData[i*nx + k] = (data[i*nx + k] - mean[i]) / std[i];
        }
    }
    for (int i = 0; i < ny; i++) {
        result[i*ny + i]=1.0;
        for (int j = 0; j < i; j++) {
            double dotProduct = 0.0;
            for (int k = 0; k < nx; k++) {
                dotProduct += normalizedData[i*nx + k] * normalizedData[j*nx + k];
            }
            result[i+j*ny] = (float)(dotProduct);
        }
    }
    delete[] mean;
    delete[] std;
    delete[] normalizedData;
}
