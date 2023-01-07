/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in in[x + y*nx]
- for each pixel (x, y), store the median of the pixels (a, b) which satisfy
  max(x-hx, 0) <= a < min(x+hx+1, nx), max(y-hy, 0) <= b < min(y+hy+1, ny)
  in out[x + y*nx].
*/
#include <vector>
#include <algorithm>

double calmean(int is, int ie, int js, int je, int nx, const float *in){
  int a=ie-is, b=je-js;
  std::vector<float> tmp(a*b);

  #pragma omp parallel
  #pragma omp for schedule(static,1) nowait
  for (int j=0; j<b; ++j){
    for (int i=0; i<a; ++i){
      tmp[i+j*a]=in[is+i+(js+j)*nx];
    }
  }
  auto n=tmp.size()/2;
  std::nth_element(tmp.begin(), tmp.begin()+n, tmp.end());
  auto med=tmp[n];
  if(!(tmp.size()&1)){
    auto med2 = std::max_element(tmp.begin(), tmp.begin()+n);
    med = (*med2 + med) / 2.0;
  }
  return med;
}

void mf(int ny, int nx, int hy, int hx, const float *in, float *out) {
  //out(ny*nx)
  #pragma omp parallel
  #pragma omp for schedule(dynamic,1) nowait
  for (int j=0; j<ny; ++j){
    for (int i=0; i<nx; ++i){
      out[i+j*nx]=calmean(std::max(i-hx, 0),std::min(i+hx+1, nx),std::max(j-hy, 0),std::min(j+hy+1, ny), nx, in);
    }
  }
}
