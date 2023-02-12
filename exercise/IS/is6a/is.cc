#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <omp.h>
struct Result
{
    int y0;
    int x0;
    int y1;
    int x1;
    float outer[3];
    float inner[3];
};
struct Result2
{
    int y0;
    int x0;
    int y1;
    int x1;
    float outer;
    float inner;
    float cost;
};
Result tor(Result2 res){
    return {res.y0, res.x0, res.y1, res.x1, {res.outer, res.outer, res.outer}, {res.inner, res.inner, res.inner}};
}
Result segment(int ny, int nx, const float *data)
{
    int nxy=nx*ny;
    std::vector<int> new_data(nxy);
    #pragma omp parallel for
    for (int y=0; y<ny; ++y){
        for (int x=0; x<nx; ++x){
            int num=x+nx*y;
            new_data[num] = int(data[3*num]);
        }
    }
    std::vector<int> new_sto(nxy+nx+ny+1,0);
    for (int y=1; y<ny+1; ++y){
        int linesum=0;
        for (int x=1; x<nx+1; ++x){
            int num=x+(nx+1)*y;
            int num2=x-1+nx*(y-1);
            linesum+=new_data[num2];
            new_sto[num]=linesum;
            new_sto[num]+=new_sto[num-nx-1];
        }
    }
    
    int all = new_sto[nxy+nx+ny];
    // Result result = {0, 0, 1, 1, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    Result2 result={0,0,1,1,0.0,0.0,INFINITY};
    // float bestcost = INFINITY;
    #pragma omp parallel
    {
        // float local_bestcost = bestcost;
        Result2 local_result = result;
#pragma omp for schedule(dynamic,1) nowait firstprivate(result)
        for (int h = 1; h <= ny; ++h)
        {
            for (int w = 1; w <= nx; ++w)
            {
                // for (int y0 = 0; y0+h <= ny; ++y0)
                for (int y1 = h; y1<=ny; ++y1)
                {
                    // for (int x0 = 0; x0+w <= nx; ++x0)
                    for (int x1=w; x1<=nx; ++x1)
                    {
                        int y0=y1-h;
                        int x0=x1-w;
                        int dim = h * w;
                        if (dim == nxy)
                            break;
                        float inv_dim = 1.0 / (dim);
                        float inv_outdim = 1.0 / (nxy - dim);
                        int n11 = x1 + (nx + 1) * y1;
                        int n01 = x0 + (nx + 1) * y1;
                        int n10 = x1 + (nx + 1) * y0;
                        int n00 = x0 + (nx + 1) * y0;
                        int inner_t = new_sto[n11] + new_sto[n00] - new_sto[n01] - new_sto[n10];
                        int outer_t = all - inner_t;
                        float inner_t_ave = (float)inner_t*inv_dim;
                        float outer_t_ave = (float)outer_t*inv_outdim;
                        float t_cost = -inner_t * inner_t_ave - outer_t * outer_t_ave;
                        if (t_cost<local_result.cost)
                        {
                            // local_result = {y0, x0, y1, x1, {outer_t_ave, outer_t_ave, outer_t_ave}, {inner_t_ave, inner_t_ave, inner_t_ave}};
                            // local_bestcost = t_cost;
                            local_result={y0,x0,y1,x1,outer_t_ave,inner_t_ave,t_cost};
                        }
                    }
                }
            }
        }
#pragma omp critical
        {
            if (local_result.cost < result.cost)
            {
                result = local_result;
                // bestcost = local_bestcost;
            }
        }
    }
    return tor(result);
}