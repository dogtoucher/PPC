#include <iostream>
#include <vector>
#include <limits>

void step(float* r, const float* d, int n){
    std::vector<float> t(n*n);
    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            t[n*j+i] = d[n*i+j];
        }
    }
    
    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            float v=std::numeric_limits<float>::infinity();
            for (int k=0; k<n; ++k){
                float x=d[n*i+k];
                float y=t[n*j+k];
                float z=x+y;
                v=std::min(v,z);
            }
            r[n*i+j]=v;
        }
    }
}