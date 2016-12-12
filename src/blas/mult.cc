
void
multiply(
  int nrows,
  double* vector,
  double* residual,
  double** matrixRows,
  int* nnzPerRow,
  int** nonzerosInRow)
{
  #pragma omp parallel for
  for (int r=0; r < nrows; ++r){
    double res = 0; 
    int nnz = nnzPerRow[r];
    int* nonzeros = nonzerosInRow[r];
    double* nonzeroVals = matrixRows[r];
    for (int i=0; i < nnz; ++i){
      int col = nonzeros[i];
      double matrixVal = nonzeroVals[i];
      double vectorVal = vector[col];
      res += matrixVal * vectorVal;
    }
    residual[r] = res;
  }
}


