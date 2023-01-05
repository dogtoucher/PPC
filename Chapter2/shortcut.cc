#include <iostream>
#include <chrono>
void step(float* r, const float* d, int n);
//generate the random number
int randInt(int min, int max){
    return min+(std::rand()%(max-min+1));
}
int main(){
    // constexpr int n;
    int n;
    std::cout << "Type the number of nodes:";
    std::cin >> n;
    // const float d[n*n] = {
    //     0, 8, 2,
    //     1, 0, 9,
    //     4, 5, 0,
    // };
    float d[n*n];
    for (int i=0; i<n*n; ++i){
        d[i]=randInt(1,9);
    }
    for (int i=0; i<n; ++i){
        d[i*(n+1)]=0;
    }
    float r[n*n];
    auto begin=std::chrono::high_resolution_clock::now();
    step(r, d, n);
    auto end=std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    // for (int i=0; i<n; ++i){
    //     for (int j=0; j<n; ++j){
    //         std::cout << r[i*n + j] << " ";
    //     }
    //     std::cout << "\n";
    // }
    std::cout << "Time taken by function: " << elapsed.count() << " nanoseconds"<< std::endl;

}