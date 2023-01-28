#include<cmath>
#include<vector>
#include<iostream>
void correlate(int ny, int nx, const float* data, float* result) {
    std::vector<double> mean;
    std::vector<double> std;
    std::vector<double> normalizedData;
    mean.reserve(ny);
    std.reserve(ny);
    normalizedData.reserve(nx*ny);

    for (int i = 0; i < ny; ++i) {
        double meani=0.0;
        for (int k = 0; k < nx; ++k) {
            meani += data[i*nx + k];
        }
        mean.push_back(meani/nx);
    }
    for (int i = 0; i < ny; ++i) {
        double sumSquares = 0.0;
        for (int k = 0; k < nx; ++k) {
            sumSquares += pow((data[i*nx + k]-mean[i]),2);
        }
        std.push_back(sqrt(sumSquares));
    }
    for (int i = 0; i < ny; ++i) {
        for (int k = 0; k < nx; ++k) {
            normalizedData.push_back((data[i*nx + k] - mean[i]) / std[i]);
        }
    }
    for (int i=0; i< ny; ++i){
        result[i*ny + i]=1.0;
    }

    for (int i = 1; i < ny; ++i) {
        for (int j = 0; j < i; ++j) {
            double dotProduct = 0.0;
            for (int k = 0; k < nx; ++k) {
                dotProduct += normalizedData[i*nx + k] * normalizedData[j*nx + k];
            }
            result[i+j*ny] = dotProduct;
        }
    }
}
