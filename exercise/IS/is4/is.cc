#include <algorithm>
#include <cmath>
#include <vector>
// #include <iostream>
typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));
// void output4(double4_t vd)
// {
//     std::cout << vd[0] << " " << vd[1] << " " << vd[2] << std::endl;
// }
constexpr double4_t d4init{
    0.0, 0.0, 0.0, 0.0};
static inline double sum4(double4_t vd)
{
    return vd[0] + vd[1] + vd[2]; //+vd[3];
}

static inline double4_t mul4i(double4_t vd, int num)
{
    double4_t tmp;
    for (int i = 0; i < 3; ++i)
    {
        tmp[i] = num * vd[i];
    }
    return tmp;
}
static inline double4_t mul4d(double4_t vd, double num)
{
    double4_t tmp;
    for (int i = 0; i < 3; ++i)
    {
        tmp[i] = num * vd[i];
    }
    return tmp;
}

struct Result
{
    int y0;
    int x0;
    int y1;
    int x1;
    float outer[3];
    float inner[3];
};

Result segment(int ny, int nx, const float *data)
{
    int nxy = nx * ny;
    // init new_data store vector
    std::vector<double4_t> new_data(nxy + nx + ny);

#pragma omp parallel for
    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
        {
            int num = x + nx * y;
            // int num2=x-1+nx*(y-1);
            new_data[num][0] = data[3 * num];
            new_data[num][1] = data[3 * num + 1];
            new_data[num][2] = data[3 * num + 2];
        }
    }
    // init summation of position xy-00
    std::vector<double4_t> new_sto(nxy + nx + ny + 1);
    for (int y = 1; y < ny + 1; ++y)
    {
        double4_t linesum = d4init;
        for (int x = 1; x < nx + 1; ++x)
        {
            int num = x + (nx + 1) * y;
            int num2 = x - 1 + nx * (y - 1);
            linesum += new_data[num2];
            new_sto[num] = linesum;
            new_sto[num] += new_sto[num - nx - 1];
        }
    }
    double4_t all = new_sto[nxy + nx + ny];

    Result result = {0, 0, 1, 1, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double bestcost = INFINITY;
// run through h and w
#pragma omp parallel shared(bestcost, result)
    {
        double local_bestcost = bestcost;
        Result local_result = result;
#pragma omp for schedule(dynamic)
        for (int h = 1; h <= ny; ++h)
        {
            for (int w = 1; w <= nx; ++w)
            {
                for (int y0 = 0; y0 <= ny - h; ++y0)
                {
                    for (int x0 = 0; x0 <= nx - w; ++x0)
                    {
                        int dim = h * w;
                        int outdim = nxy - dim;
                        if (dim == nxy)
                            break;
                        double inv_dim = 1.0 / (dim);
                        double inv_outdim = 1.0 / (outdim);
                        int x1 = x0 + w;
                        int y1 = y0 + h;
                        int n11 = x1 + (nx + 1) * y1;
                        int n01 = x0 + (nx + 1) * y1;
                        int n10 = x1 + (nx + 1) * y0;
                        int n00 = x0 + (nx + 1) * y0;
                        double4_t inner_t = new_sto[n11] + new_sto[n00] - new_sto[n01] - new_sto[n10];
                        double4_t outer_t = all - inner_t;
                        // inner and outer[c]
                        auto inner_t_ave = mul4d(inner_t, inv_dim);
                        auto outer_t_ave = mul4d(outer_t, inv_outdim);
                        double4_t t_cost4 = -inner_t * inner_t_ave - outer_t * outer_t_ave;
                        double t_cost = sum4(t_cost4);
                        if (t_cost < local_bestcost)
                        {
                            local_result = {y0, x0, y1, x1, {(float)outer_t_ave[0], (float)outer_t_ave[1], (float)outer_t_ave[2]}, {(float)inner_t_ave[0], (float)inner_t_ave[1], (float)inner_t_ave[2]}};
                            local_bestcost = t_cost;
                        }
                    }
                }
            }
        }
#pragma omp critical
        {
            if (local_bestcost < bestcost)
            {
                result = local_result;
                bestcost = local_bestcost;
            }
        }
    }
    return result;
}