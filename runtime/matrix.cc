#include <random>
#include <mkl.h>

#include "matrix.h"

using namespace std;

void Matrix::randomFill(){
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0,1);
  int nrows = blockSize_ * blockGridSize_;
#pragma omp parallel for
  for(int i = 0; i < nrows * nrows; i++){
    //double x = double(0.05 + (i*7)%4 - (i*5)%2) / ((i*11)%5+0.1);
    double x = dis(gen);
    storage_[i] = x;
  }
}

void Matrix::symmetricFill(){
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0,1);
  int nrows = blockSize_ * blockGridSize_;
#pragma omp parallel for
  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < nrows; j++){
      int idx = i*nrows + j;
      if(i == j){
        double x = dis(gen);
        storage_[idx] = x;
      } else {
        storage_[idx] = 0.0;
      }
    }
  }
}

Matrix::Matrix(int matrixSize, int blockSize) : 
  blockSize_(blockSize), 
  blockGridSize_(matrixSize)
{
  storage_ = (double*) mkl_malloc(matrixSize * matrixSize * blockSize * blockSize * sizeof(double), 64);
}

void Matrix::print(const char* label){
  if (label) printf("%s\n", label);
  int nrows = blockSize_ * blockGridSize_;
  for(int i = 0; i < nrows; i++){
    if(i > 0 && i % blockSize_ == 0){
      printf("\n");
    }
    for(int j = 0; j < nrows; j++){
      if(j > 0 && j % blockSize_ == 0){
        printf("   ");
      }
      auto val = storage_[j*nrows + i];
      printf("%6.4f ", val);
    }
    printf("\n");
  }
  printf("\n");

}
