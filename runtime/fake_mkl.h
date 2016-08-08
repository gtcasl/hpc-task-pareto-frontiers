#include <cblas.h>


extern "C"
void
vdSub(int, double*, double*, double*);

extern "C"
void 
cblas_daxpby (const int n, const double a, const double *x, 
  const int incx, const double b, double *y, const int incy);

extern "C"
void
dgetrf(int*, int*, double*, int*, int*, int*);

extern "C"
void
dpotrf(char*, int*, double*, int*, int*);

inline int
myDPOTRF(char ty, int N, double* A, int LDA){
  int info;
  dpotrf(&ty, &N, A, &LDA, &info);
  return info;
}

#define LAPACKE_dpotrf(ignore, ty, size1, A, size2) \
  myDPOTRF(ty, size1, A, size2)

