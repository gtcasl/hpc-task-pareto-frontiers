#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
#include <unistd.h>
#include <vector>

#include "mkl.h"
#include <omp.h>

#include <lwperf/lwperf.h>

// compile with -mkl
// link with -lmkl_rt

using namespace std;

void rand_fill(double* a,int n){
  for(int i = 0; i < n; i++){
    a[i] = ((double)rand()/(double)RAND_MAX);
  }
}

struct mklfreer{
  void operator()(double* d){
    mkl_free((void*) d);
  }
};

struct mklvec{
  unique_ptr<double,mklfreer> v;
  mklvec(int n){
    v.reset((double*)mkl_malloc(n * sizeof(double), 64));
  }
  double* data(){
    return v.get();
  }
};

// Cholesky
// DPOTRF DTRSM DSYRK DGEMM
struct dpotrf_t{
  mklvec a;
  int n;
  dpotrf_t(int N) : a(N*N), n(N){
    rand_fill(a.data(), N*N);
    mklvec a_copy(N*N);
    for(int i = 0; i < N*N; i++){
      a_copy.data()[i] = a.data()[i];
    }
    mklvec sym_mat(N*N);
    mkl_domatadd('r', 'n', 't', N, N, 1.0, a.data(), N, 1.0, a_copy.data(), N,
                 sym_mat.data(), N);

    mklvec ident(N*N);
    for(int i = 0; i < N*N; i++){
      ident.data()[i] = 0.0;
    }
    for(int i = 0; i < N; i++){
      ident.data()[i * N + i] = N;
    }
    mkl_domatadd('r', 'n', 'n', N, N, 1.0, a_copy.data(), N, 1.0, ident.data(),
                 N, a.data(), N);
  }
};

void run_dpotrf(dpotrf_t& d){
  mklvec a(d.n*d.n);
  cblas_dcopy(d.n*d.n, d.a.data(), 1, a.data(), 1);
  int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', d.n, a.data(), d.n);
  if(info != 0){
    cerr << "Error running dpotrf" << endl;
    exit(-1);
  }
}

struct dtrsm_t{
  mklvec a,b;
  int n;
  dtrsm_t(int N) : a(N*N), b(N*N), n(N){
    rand_fill(a.data(), N*N);
    rand_fill(b.data(), N*N);
  }
};

void run_dtrsm(dtrsm_t& d){
  cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasUnit,
              d.n, d.n, 1.0, d.a.data(), d.n, d.b.data(), d.n);
}

struct dsyrk_t{
  mklvec a,c;
  int n;
  dsyrk_t(int N) : a(N*N), c(N*N), n(N){
    rand_fill(a.data(), N*N);
    rand_fill(c.data(), N*N);
  }
};

void run_dsyrk(dsyrk_t& d){
  cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, d.n, d.n, 1.0,
              d.a.data(), d.n, 1.0, d.c.data(), d.n);
}

struct dgemm_t{
  mklvec a,b,c;
  int n;
  dgemm_t(int N) : a(N*N), b(N*N), c(N*N), n(N){
    rand_fill(a.data(), N*N);
    rand_fill(b.data(), N*N);
  }
};

void run_dgemm(dgemm_t& d){
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, d.n, d.n, d.n, 1.0,
              d.a.data(), d.n, d.b.data(), d.n, 0, d.c.data(), d.n);
}

// CG
// DGEMV DAXPY DDOT
struct dgemv_t{
  mklvec x,a,y;
  int n;
  dgemv_t(int N) : x(N), a(N*N), y(N), n(N){
    // initialize data with garbage or whatever
    rand_fill(x.data(), N);
    rand_fill(a.data(), N*N);
  }
};

void run_dgemv(dgemv_t& d){
  cblas_dgemv(CblasRowMajor, CblasNoTrans, d.n, d.n, 1, d.a.data(), d.n,
              d.x.data(), 1, 0, d.y.data(), 1);
}

struct daxpy_t{
  mklvec x,y;
  double a;
  int n;
  daxpy_t(int N) : x(N*N), y(N*N), a(43.59), n(N*N){
    rand_fill(x.data(), N*N);
    rand_fill(y.data(), N*N);
  }
};

void run_daxpy(daxpy_t& d){
  cblas_daxpy(d.n, d.a, d.x.data(), 1, d.y.data(), 1);
}

struct ddot_t{
  mklvec x,y;
  int n;
  ddot_t(int N) : x(N*N), y(N*N), n(N*N){
    rand_fill(x.data(), N*N);
    rand_fill(y.data(), N*N);
  }
};

void run_ddot(ddot_t& d){
  cblas_ddot(d.n, d.x.data(), 1, d.y.data(), 1);
}

template<typename T>
void experiment(lwperf_t perf, const char* name, void(*run_fn)(T&), int N,
                int num_iters){
  cout << "Running experiment for " << name << endl;
  T prepped(N);

  // warm up
  run_fn(prepped);
  lwperf_log(perf, name);
  for(int i = 0; i < num_iters; i++){
    run_fn(prepped);
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
#define EXPERIMENT(name,iters) experiment(perf, #name, run_ ## name, N, iters)
  EXPERIMENT(dpotrf,40);
  EXPERIMENT(dtrsm,800);
  EXPERIMENT(dsyrk,160);
  EXPERIMENT(dgemm,400);
  EXPERIMENT(dgemv,8000);
  EXPERIMENT(daxpy,40000);
  EXPERIMENT(ddot,40000);
#undef EXPERIMENT

  lwperf_finalize(perf);
  return 0;
}

