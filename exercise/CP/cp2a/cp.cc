/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/

#include<cmath>
#include <vector>
void correlate(int ny, int nx, const float *data_, float *result) {
// cal the corr between row i and row j
    constexpr int nb = 8;
    int na = (nx-1+nb)/nb;
    int nab = na*nb;
    std::vector<double> data(ny*nab, 0.0);

    for (int i = 0; i < ny; ++i) {
        for (int k = 0; k < nx; ++k) {
            data[k + i*nab]=data_[k + i*nx];
        }
    }   

    std::vector<double> total(ny, 0.0);
    std::vector<double> total2(ny, 0.0);
    std::vector<double> bottom(ny, 0.0);
    for (int i=0; i<ny; ++i){
        result[i + i*ny] = 1;
        for (int k=0; k<nx; ++k){
            total[i] += double(data[k + i*nab]);
            total2[i] += double(data[k + i*nab])*double(data[k + i*nab]);
        }
        bottom[i] = sqrt(total2[i]*nx-total[i]*total[i]);
    }

    for (int i=0; i<ny; ++i){
        for (int j=0; j<i; ++j){
            double total_ij[nb]={0.0};
            for (int ka=0; ka<na; ++ka){
                for (int kb=0; kb<nb; ++kb){
                    total_ij[kb] += double(data[kb + nb*ka + i*nab])*double(data[kb + nb*ka + j*nab]);
                }   
            }

            double total_ij_ = 0.0;
            for (int k=0; k<nb; k++){
                total_ij_+=total_ij[k];
            }
            double top=total_ij_*nx-(total[i]*total[j]);
            double bottom_=bottom[i]*bottom[j];
            result[i + j*ny]=float(top/bottom_);
        }

    }
}

