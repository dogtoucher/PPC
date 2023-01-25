#include <algorithm>
#include <omp.h>
#include <cmath>

typedef unsigned long long data_t;

void merge_abc(int a, int b, int c, data_t *data){
    std::inplace_merge(data+a, data+b, data+c);
}
void sort_ab(int a, int b, int layer, data_t *data, int num, int n){
    if(layer==0){
        std::sort(data+a, data+b);
    }
    else{
        int tmp=(b-a)/num;
        for(int i=0; i<num; ++i){
            // sort_ab(a+i*tmp, std::min(a+(i+1)*tmp,b), layer-1, data, num, n);
            //i*num ~ min(n,(i+1)*num) sorted
        }
        std::sort(data+a, data+b);
    }
}
void psort(int n, data_t *data){
    // std::sort(data, data+n);
    int num=omp_get_max_threads();
    // int num=8;
    int layer=log(n)/log(num);
    sort_ab(0,n,layer,data,num,n);
}