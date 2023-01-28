#include <omp.h>
#include <cmath>
#include <algorithm>
typedef unsigned long long data_t;

void ppsort(int t, int n, data_t *data) {
    if(t<=0){
        std::sort(data, data+n);
        return;
    }
    int m = (n+1)/2;
    #pragma omp task
    ppsort(t-1, m, data);
    #pragma omp task
    ppsort(t-1, n-m, data+m);
    #pragma omp taskwait
    std::inplace_merge(data, data+m, data+n);
}
void psort(int n, data_t *data) {
    int t=std::log2(omp_get_max_threads())*2;
    #pragma omp parallel
    {
        #pragma omp single
        {
            ppsort(t, n, data);
        }
    }
}


