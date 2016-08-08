#ifdef no_mkl

extern "C" {

void
vdSub(int n, double* a, double* b, double* c){
#pragma omp parallel for
  for (int i=0; i < n; ++i){
    c[i] = a[i] - b[i];
  }
}

void 
cblas_daxpby (const int n, const double a, const double *x, 
  const int incx, const double b, double *y, const int incy){
#pragma omp parallel for
  for (int i=0; i < n; ++i){
    y[i] = a*x[i] + b*y[i];
  }
}

}


#endif

