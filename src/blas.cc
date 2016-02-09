#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <unistd.h>
#include <vector>

#include "mkl.h"
#include <omp.h>

#include <lwperf/lwperf.h>

// compile with -mkl
// link with -lmkl_rt

#define NUM_ITERS 50

using namespace std;

void rand_fill(double* a,int n){
  for(int i = 0; i < n; i++){
    a[i] = ((double)rand()/(double)RAND_MAX);
  }
}

double gen_rand(){
  return ((double)rand()/(double)RAND_MAX);
}

// Cholesky
// DPOTRF DTRSM DSYRK DGEMM
struct dpotrf_t{
  vector<double> a;
  dpotrf_t(int N) : a(N * N){
    generate(a.begin(), a.end(), gen_rand);
    vector<double> a_copy = a;
    vector<double> sym_mat(N * N);
    mkl_domatadd('r', 'n', 't', N, N, 1.0, &a.front(), N, 1.0,
                 &a_copy.front(), N, &sym_mat.front(), N);

    vector<double> ident(N*N, 0.0);
    for(int i = 0; i < N; i++){
      ident[i * N + i] = N;
    }
    mkl_domatadd('r', 'n', 'n', N, N, 1.0, &a_copy.front(), N, 1.0,
                 &ident.front(), N, &a.front(), N);
  }
};

void run_dpotrf(dpotrf_t d){
  int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', d.a.size(), &d.a.front(),
                            d.a.size());
  if(info != 0){
    cerr << "Error running dpotrf" << endl;
    exit(-1);
  }
}

struct dtrsm_t{
  vector<double> a,b;
  dtrsm_t(int N) : a(N*N), b(N*N){
    generate(a.begin(), a.end(), gen_rand);
    generate(b.begin(), b.end(), gen_rand);
  }
};

void run_dtrsm(dtrsm_t d){
  cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasUnit,
              d.a.size(), d.a.size(), 1.0, &d.a.front(), d.a.size(),
              &d.b.front(), d.a.size());
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

template<typename T>
void experiment(lwperf_t perf, const char* name, void(*run_fn)(T), int N){
  cout << "Running experiment for " << name << endl;
  vector<T> prepped;
  for(int i = 0; i < NUM_ITERS + 1; i++){
    prepped.push_back(T(N));
  }

  // warm up
  run_fn(prepped[NUM_ITERS]);
  lwperf_log(perf, name);
  for(int i = 0; i < NUM_ITERS; i++){
    run_fn(prepped[i]);
  }
  lwperf_stop(perf, name);
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

// DPOTRF DTRSM DSYRK DGEMM
// DGEMV DAXPY DDOT
#define EXPERIMENT(name) experiment(perf, #name, run_ ## name, N)
  EXPERIMENT(dpotrf);
#undef EXPERIMENT

  lwperf_finalize(perf);
  return 0;
}

