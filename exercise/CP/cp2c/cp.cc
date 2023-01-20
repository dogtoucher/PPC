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

typedef double double4_t __attribute__ ((vector_size (4*sizeof(double))));
void correlate(int ny, int nx, const float *data_, float *result) {
// cal the corr between row i and row j
    constexpr int nb = 4;
    int na = (nx-1+nb)/nb;

    std::vector<double4_t> vdata(ny*na);

    for (int i = 0; i < ny; ++i) {
        for (int ka = 0; ka < na; ++ka) {
            for (int kb = 0; kb < nb; ++kb){
                int j = ka*nb+kb;
                vdata[na*i+ka][kb] = j<nx ? data_[j+i*nx] : 0.0;
            }
        }
    }   

    std::vector<double> total(ny, 0.0);
    std::vector<double> total2(ny, 0.0);
    std::vector<double> bottom(ny, 0.0);

    for (int i=0; i<ny; ++i){
        result[i + i*ny] = 1;
    }
    
    for (int i=0; i<ny; ++i){
        for (int ka=0; ka<na; ++ka){
            for (int kb=0; kb<nb; ++kb){
                total[i] += vdata[na*i+ka][kb];
                total2[i] += vdata[na*i+ka][kb]*vdata[na*i+ka][kb];
            }
        }
        bottom[i] = sqrt(total2[i]*nx-total[i]*total[i]);
    }

    for (int i=0; i<ny; ++i){
        for (int j=0; j<i; ++j){
            // std::vector<double> total_ij(nb, 0.0);
            double total_ij[nb];//={0.0};
            for (int kb=0; kb<nb; ++kb){
                total_ij[kb]=0.0;
            }
            for (int ka=0; ka<na; ++ka){
                for (int kb=0; kb<nb; ++kb){
                    total_ij[kb] += vdata[na*i+ka][kb]*vdata[na*j+ka][kb];
                }   
            }
            double total_ij_ = 0.0;
            for (int kb=0; kb<nb; ++kb){
                total_ij_+=total_ij[kb];
            }
            double top=total_ij_*nx-(total[i]*total[j]);
            double bottom_=bottom[i]*bottom[j];
            result[i + j*ny]=float(top/bottom_);
        }
    }
}

