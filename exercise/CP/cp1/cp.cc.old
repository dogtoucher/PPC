/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
#include<cmath>
#include<iostream>
double corr(int i, int j, int nx, const float *data){
    // if (i==j) return 1.0;
    double total_i=0, total_j=0;
    double total_ii=0, total_jj=0;
    double total_ij=0;
    for (int k=0; k<nx; k++){
        total_i += double(data[k + i*nx]);
        total_j += double(data[k + j*nx]);

        total_ii += double(data[k + i*nx])*double(data[k + i*nx]);
        total_jj += double(data[k + j*nx])*double(data[k + j*nx]);

        total_ij += double(data[k + i*nx])*double(data[k + j*nx]);
    }
    double top=total_ij-(total_i*total_j/nx);
    double bottom=sqrt((total_ii-pow(total_i,2)/nx)*(total_jj-pow(total_j,2)/nx));
    if (bottom==0) return 0;
    return top/bottom;
}


void correlate(int ny, int nx, const float *data, float *result) {
// cal the corr between row i and row j
    for (int i=0; i<ny; ++i){
        result[i + i*ny] = 1;
        for (int j=0; j<i; ++j){
            result[i + j*ny] = float(corr(i,j,nx,data));
        }

    }

}
