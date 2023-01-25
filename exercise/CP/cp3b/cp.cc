#include<cmath>
#include<vector>

typedef float float8_t __attribute__ ((vector_size (8*sizeof(float))));
constexpr float8_t f8init {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

static inline float sum8(float8_t vd) {
    return vd[0]+vd[1]+vd[2]+vd[3]+vd[4]+vd[5]+vd[6]+vd[7];
}

void correlate(int ny, int nx, const float* data, float* result) {
    constexpr int nb = 8;
    int na = (nx-1+nb)/nb;

    constexpr int nd = 8;
    int nc = (ny+nd-1)/nd;
    int ncd = nc*nd;
    std::vector<float8_t> vdata(ncd*na);

    std::vector<float> mean(ny,0.0);
    std::vector<float> std(ny,0.0);

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < ny; ++i) {
        for (int k = 0; k < nx; ++k) {
            mean[i] += data[i*nx+k];
        }
        mean[i] /= nx;
    }
    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < ny; ++i) {
        float sumSquares = 0.0;
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
            vdata[na*i+ka] = f8init;
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
            float8_t tmp[nd][nd];
            float8_t x[nd];
            float8_t y[nd];
            for(int id = 0; id < nd; ++id){
                for(int jd = 0; jd < nd; ++jd){
                    tmp[id][jd] = f8init;
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
            result[j+i*ny] = sum8(tmp[id][jd]);
          }
        }
      }
    }
  }
}