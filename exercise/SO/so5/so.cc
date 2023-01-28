#include <algorithm>
#include <omp.h>
#include <cmath>
typedef unsigned long long data_t;

void ppsort(int t, int n, data_t *data) {

    if (n < 2) return;
    if(t<=0){
        std::sort(data, data+n);
        return;
    }
    auto p1 = data[n/2];
    auto p2 = data[0];
    auto p3 = data[n-1];
    auto pivot = p1+p2+p3-std::max(p1, std::max(p2,p3))-std::min(p1, std::max(p2,p3));
    auto *left = data;
    auto *right = data + n - 1;

    while (left <= right) {
        while (*left < pivot) left++;
        while (*right > pivot) right--;
        if (left <= right) {
            auto temp = *left;
            *left = *right;
            *right = temp;
            left++;
            right--;
        }
    }
    #pragma omp task
    ppsort(t-1, right - data + 1, data);
    #pragma omp task
    ppsort(t-1, n - (left - data), left);
}

void psort(int n, data_t *data) {
    int t=int(std::log2(omp_get_max_threads()))*2;
    #pragma omp parallel
    {
        #pragma omp single
        {
            ppsort(t, n, data);
        }
    }
}
