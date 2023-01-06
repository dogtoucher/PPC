/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/

#include<cmath>

void correlate(int ny, int nx, const float *data, float *result) {
// cal the corr between row i and row j
    double total[ny]={0.0};
    double total2[ny]={0.0};
    double bottom[ny]={0.0};
    for (int i=0; i<ny; ++i){
        result[i + i*ny] = 1;
        for (int k=0; k<nx; ++k){
            total[i] += double(data[k + i*nx]);
            total2[i] += double(data[k + i*nx])*double(data[k + i*nx]);
        }
        bottom[i] = sqrt(total2[i]*nx-total[i]*total[i]);
    }
    for (int i=0; i<ny; ++i){
        for (int j=0; j<i; ++j){
            double total_ij=0.0;
            for (int k=0; k<nx; ++k){
                total_ij += double(data[k + i*nx])*double(data[k + j*nx]);
            }
            double top=total_ij*nx-(total[i]*total[j]);
            double bottom_=bottom[i]*bottom[j];
            result[i + j*ny]=float(top/bottom_);
        }

    }
}

