#include<cmath>
#include<vector>

typedef double double4_t __attribute__ ((vector_size (4*sizeof(double))));
constexpr double4_t d4init {
    0.0, 0.0, 0.0, 0.0
};

static inline double sum4(double4_t vd) {
    return vd[0]+vd[1]+vd[2]+vd[3];
}

void correlate(int ny, int nx, const float* data, float* result) {
    constexpr int nb = 4;
    int na = (nx-1+nb)/nb;

    constexpr int nd = 8;
    int nc = (ny+nd-1)/nd;
    int ncd = nc*nd;
    std::vector<double4_t> vdata(ncd*na);

    std::vector<double> mean(ny,0.0);
    std::vector<double> std(ny,0.0);

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < ny; ++i) {
        for (int k = 0; k < nx; ++k) {
            mean[i] += data[i*nx+k];
        }
        mean[i] /= nx;
    }
    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < ny; ++i) {
        double sumSquares = 0.0;
        for (int k = 0; k < nx; ++k) {
            sumSquares += pow(data[i*nx+k]-mean[i],2);
        }
        std[i] = sqrt(sumSquares);
    }

    // #pragma omp parallel for schedule(static,1)
    // for (int i = 0; i < ny; ++i) {
    //     mean[i] /= std[i];
    // }

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < ny; ++i) {
        for (int ka = 0; ka < na; ++ka){
            for (int kb=0; kb < nb; ++kb){
                int j = ka*nb+kb;
                vdata[na*i+ka][kb] = j<nx ? (data[i*nx+j] - mean[i])/std[i] : 0.0;
            }
        }
    }

    #pragma omp parallel for schedule(static,1)
    for (int i = ny; i < ncd; ++i){
        for (int ka = 0; ka < na; ++ka){
            vdata[na*i+ka] = d4init;
        }
    }
    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < ny; ++i) {
        result[i*ny+i]=1.0;
    }

    #pragma omp parallel 
    #pragma omp for schedule(dynamic,1) nowait
    for(int ic = 0; ic < nc; ++ic){
        for(int jc = ic; jc < nc; ++jc){
            double4_t tmp[nd][nd];
            double4_t x[nd];
            double4_t y[nd];
            for(int id = 0; id < nd; ++id){
                for(int jd = 0; jd < nd; ++jd){
                    tmp[id][jd] = d4init;
                }
            } 
      for(int ka = 0; ka < na; ++ka){
        for (int kd=0; kd < nd; ++kd){
            x[kd]=vdata[na*(ic*nd+kd)+ka];
            y[kd]=vdata[na*(jc*nd+kd)+ka];
        }
        for(int id = 0; id < nd; ++id){
            for(int jd = 0; jd < nd; ++jd){
                tmp[id][jd] += x[id]*y[jd];
            }
        } 
      }
      for(int id = 0; id < nd; ++id){
        for(int jd = 0; jd < nd; ++jd){
          int i = ic*nd+id;
          int j = jc*nd+jd;
          if(i<j && i<ny && j<ny){
            result[j+i*ny] = (float)(sum4(tmp[id][jd]));
          }
        }
      }
    }
  }
}