#include <cg/types.h>

void
multiply(
  int nrows,
  DoubleArray vector,
  DoubleArray residual,
  DoubleArray matrix,
  IntArray nnzPerRow,
  IntArray nonzerosInRow)
{
  double* v = vector;
  double* r = residual;
  double* A = matrix;
  int* nonzeros = nonzerosInRow;
  #pragma omp parallel for
  for (int row=0; row < nrows; ++row){
    double res = 0; 
    int nnz = nnzPerRow[row];
    for (int i=0; i < nnz; ++i){
      int col = nonzeros[i];
      res += A[i] * v[col];
    }
    A += nnz;
    nonzeros += nnz;
    r[row] = res;
  }
}


