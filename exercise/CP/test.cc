#include <iostream>
#include <cmath>

using namespace std;

// Function to calculate Pearson correlation coefficient between rows i and j
double pearson_correlation(int i, int j, int nx, double data[2][2])
{
    double sum_i = 0, sum_j = 0;
    double sum_sq_i = 0, sum_sq_j = 0;
    double p_sum = 0;

    for (int k = 0; k < nx; k++)
    {
        sum_i += data[i][k];
        sum_j += data[j][k];

        sum_sq_i += data[i][k] * data[i][k];
        sum_sq_j += data[j][k] * data[j][k];

        p_sum += data[i][k] * data[j][k];
    }
    cout << "sum_i" << sum_i << endl;
    cout << "sum_j" << sum_j << endl;
    cout << "p_sum" << p_sum << endl;

    double num = p_sum - (sum_i * sum_j / nx);
    double den = sqrt((sum_sq_i - pow(sum_i, 2) / nx) * (sum_sq_j - pow(sum_j, 2) / nx));

    if (den == 0) return 0;

    return num / den;
}

int main()
{
    // Example data
    // double data[5][5] = {
    //     {1, 2, 3, 4, 5},
    //     {2, 3, 4, 5, 6},
    //     {3, 4, 5, 6, 7},
    //     {4, 5, 6, 7, 8},
    //     {5, 6, 7, 8, 9}
    // };
    double data[2][2]={{+0.81472367, +0.90579194},{+0.45150527, +0.49610928}};   
    // Calculate Pearson correlation coefficient between rows 0 and 1
    double r = pearson_correlation(0, 1, 2, data);

    cout << "Pearson correlation coefficient: " << float(r) << endl;

    return 0;
}