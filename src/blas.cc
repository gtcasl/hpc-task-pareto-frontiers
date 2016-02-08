#include <cstdlib>
#include <ctime>
#include <iostream>
#include <unistd.h>
#include "mkl.h"
#include <lwperf/lwperf.h>
#include <omp.h>

// compile with -mkl
// link with -lmkl_rt

using namespace std;

void rand_fill(double* a,int n){
  for(int i = 0; i < n; i++){
    a[i] = ((double)rand()/(double)RAND_MAX);
  }
}

// Cholesky
// DPOTRF DTRSM DSYRK DGEMM
struct dpotrf_t{
  int n;
  double* a;
};

dpotrf_t prep_dpotrf(int N){
  dpotrf_t d;
  d.n = N;
  d.a = new double[d.n * d.n];
  rand_fill(d.a, d.n * d.n);
  double* a_copy = new double[d.n * d.n];
  for(int i = 0; i < d.n * d.n; i++){
    a_copy[i] = d.a[i];
  }
  double* sym_mat = new double[d.n * d.n];
  mkl_domatadd('r', 'n', 't', d.n, d.n, 1.0, d.a, d.n, 1.0, a_copy, d.n, sym_mat, d.n);

  double* ident = new double[d.n * d.n];
  for(int i = 0; i < d.n * d.n; i++){
    ident[i] = 0;
  }
  for(int i = 0; i < d.n; i++){
    ident[i * d.n + i] = d.n;
  }
  mkl_domatadd('r', 'n', 'n', d.n, d.n, 1.0, a_copy, d.n, 1.0, ident, d.n, d.a, d.n);
  delete[] a_copy;
  delete[] ident;
  return d;
}
void run_dpotrf(dpotrf_t d){
  int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', d.n, d.a, d.n);
  if(info != 0){
    cerr << "Error running dpotrf" << endl;
    exit(-1);
  }
}

struct dtrsm_t{
  int n;
  double* a;
  double* b;
};

dtrsm_t prep_dtrsm(int N){
  dtrsm_t d;
  d.n = N;
  d.a = new double[d.n * d.n];
  d.b = new double[d.n * d.n];
  rand_fill(d.a, d.n * d.n);
  rand_fill(d.b, d.n * d.n);
  return d;
}
void run_dtrsm(dtrsm_t d){
  cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasUnit, d.n, d.n, 1.0, d.a, d.n, d.b, d.n);
}

struct dsyrk_t{
  int n;
  int k;
  double* a;
  double* c;
};

dsyrk_t prep_dsyrk(int N){
  dsyrk_t d;
  d.n = N;
  d.a = new double[d.n * d.n];
  d.c = new double[d.n * d.n];
  rand_fill(d.a, d.n * d.n);
  rand_fill(d.c, d.n * d.n);
  return d;
}
void run_dsyrk(dsyrk_t d){
  cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, d.n, d.n, 1.0, d.a, d.n, 1.0, d.c, d.n);
}

struct dgemm_t{
  int n;
  double* a;
  double* b;
  double* c;
};

dgemm_t prep_dgemm(int N){
  dgemm_t d;
  d.n = N;
  d.a = new double[d.n * d.n];
  d.b = new double[d.n * d.n];
  d.c = new double[d.n * d.n];
  rand_fill(d.a, d.n * d.n);
  rand_fill(d.b, d.n * d.n);
  return d;
}
void run_dgemm(dgemm_t d){
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, d.n, d.n, d.n, 1.0, d.a, d.n, d.b, d.n, 0, d.c, d.n);
}

// CG
// DGEMV DAXPY DDOT
struct dgemv_t{
  double* x;
  double* a;
  double* y;
  int n;
};

dgemv_t prep_dgemv(int N){
  dgemv_t d;
  d.n = N;
  d.x = new double[d.n];
  d.a = new double[d.n*d.n];
  d.y = new double[d.n];
  // initialize data with garbage or whatever
  rand_fill(d.x,d.n);
  rand_fill(d.a,d.n*d.n);
  return d;
}
void run_dgemv(dgemv_t d){
  //cblas_dgemv(Layout, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
  cblas_dgemv(CblasRowMajor, CblasNoTrans, d.n, d.n, 1, d.a, d.n, d.x, 1, 0, d.y, 1);
}

struct daxpy_t{
  double *x;
  double *y;
  double a;
  int n;
};

daxpy_t prep_daxpy(int N){
  daxpy_t d;
  d.n = N * N;
  d.x = new double[d.n];
  d.y = new double[d.n];
  d.a = 43.59;
  // initialize arrays/pick a
  rand_fill(d.x,d.n);
  rand_fill(d.y,d.n);
  return d;
}
void run_daxpy(daxpy_t d){
  cblas_daxpy(d.n, d.a, d.x, 1, d.y, 1);
}

struct ddot_t{
  double *x;
  double *y;
  int n;
};

ddot_t prep_ddot(int N){
  ddot_t d;
  d.n = N * N;
  d.x = new double[d.n];
  d.y = new double[d.n];
  rand_fill(d.x,d.n);
  rand_fill(d.y,d.n);
  return d;
}
void run_ddot(ddot_t d){
  cblas_ddot(d.n, d.x, 1, d.y, 1);
}

int main(int argc, char* argv[]){
  srand((unsigned)time(NULL));

  if(argc != 2){
    cerr << "Error: Please provide size N" << endl;
    return -1;
  }
  int N = atoi(argv[1]);

  lwperf_t perf = lwperf_init("phi", "test", "database");
  lwperf_init_papi(perf);
  int nthreads = omp_get_max_threads();
  std::cout << "nthreads: " << nthreads << std::endl;
  lwperf_add_invariant(perf, "nthreads", nthreads);

  // Prep Cholesky
  cout << "Prepping dpotrf" << endl;
  dpotrf_t potrft = prep_dpotrf(N);
  cout << "Prepping dtrsm" << endl;
  dtrsm_t trsmt = prep_dtrsm(N);
  cout << "Prepping dsyrk" << endl;
  dsyrk_t syrkt = prep_dsyrk(N);
  cout << "Prepping dgemm" << endl;
  dgemm_t gemmt = prep_dgemm(N);

  // Prep CG
  cout << "Prepping dgemv" << endl;
  dgemv_t gemvt = prep_dgemv(N);
  cout << "Prepping daxpy" << endl;
  daxpy_t axpyt = prep_daxpy(N);
  cout << "Prepping ddot" << endl;
  ddot_t dott = prep_ddot(N);

  // Cholesky
  cout << "Running dpotrf" << endl;
  lwperf_log(perf, "dpotrf");
  run_dpotrf(potrft);
  lwperf_stop(perf, "dpotrf");

  cout << "Running dtrsm" << endl;
  lwperf_log(perf, "dtrsm");
  run_dtrsm(trsmt);
  lwperf_stop(perf, "dtrsm");

  cout << "Running dsyrk" << endl;
  lwperf_log(perf, "dsyrk");
  run_dsyrk(syrkt);
  lwperf_stop(perf, "dsyrk");

  cout << "Running dgemm" << endl;
  lwperf_log(perf, "dgemm");
  run_dgemm(gemmt);
  lwperf_stop(perf, "dgemm");

  // CG
  cout << "Running dgemv" << endl;
  lwperf_log(perf, "dgemv");
  run_dgemv(gemvt);
  lwperf_stop(perf, "dgemv");

  cout << "Running daxpy" << endl;
  lwperf_log(perf, "daxpy");
  run_daxpy(axpyt);
  lwperf_stop(perf, "daxpy");

  cout << "Running ddot" << endl;
  lwperf_log(perf, "ddot");
  run_ddot(dott);
  lwperf_stop(perf, "ddot");

  lwperf_finalize(perf);
  return 0;
}

