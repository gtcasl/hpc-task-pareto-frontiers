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
  int n;
  dpotrf_t(int N) : a(N * N), n(N){
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
  int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', d.n, &d.a.front(), d.n);
  if(info != 0){
    cerr << "Error running dpotrf" << endl;
    exit(-1);
  }
}

struct dtrsm_t{
  vector<double> a,b;
  int n;
  dtrsm_t(int N) : a(N*N), b(N*N), n(N){
    generate(a.begin(), a.end(), gen_rand);
    generate(b.begin(), b.end(), gen_rand);
  }
};

void run_dtrsm(dtrsm_t d){
  cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasUnit,
              d.n, d.n, 1.0, &d.a.front(), d.n, &d.b.front(), d.n);
}

struct dsyrk_t{
  vector<double> a,c;
  int n;
  dsyrk_t(int N) : a(N*N), c(N*N), n(N){
    generate(a.begin(), a.end(), gen_rand);
    generate(c.begin(), c.end(), gen_rand);
  }
};

void run_dsyrk(dsyrk_t d){
  cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, d.n, d.n, 1.0,
              &d.a.front(), d.n, 1.0, &d.c.front(), d.n);
}

struct dgemm_t{
  vector<double> a,b,c;
  int n;
  dgemm_t(int N) : a(N*N), b(N*N), c(N*N), n(N){
    generate(a.begin(), a.end(), gen_rand);
    generate(b.begin(), b.end(), gen_rand);
  }
};

void run_dgemm(dgemm_t d){
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, d.n, d.n, d.n, 1.0,
              &d.a.front(), d.n, &d.b.front(), d.n, 0, &d.c.front(), d.n);
}

// CG
// DGEMV DAXPY DDOT
struct dgemv_t{
  vector<double> x,a,y;
  int n;
  dgemv_t(int N) : x(N), a(N*N), y(N), n(N){
    // initialize data with garbage or whatever
    generate(x.begin(), x.end(), gen_rand);
    generate(a.begin(), a.end(), gen_rand);
  }
};

void run_dgemv(dgemv_t d){
  cblas_dgemv(CblasRowMajor, CblasNoTrans, d.n, d.n, 1, &d.a.front(), d.n,
              &d.x.front(), 1, 0, &d.y.front(), 1);
}

struct daxpy_t{
  vector<double> x,y;
  double a;
  int n;
  daxpy_t(int N) : x(N*N), y(N*N), a(43.59), n(N){
    generate(x.begin(), x.end(), gen_rand);
    generate(y.begin(), y.end(), gen_rand);
  }
};

void run_daxpy(daxpy_t d){
  cblas_daxpy(d.n, d.a, &d.x.front(), 1, &d.y.front(), 1);
}

struct ddot_t{
  vector<double> x,y;
  int n;
  ddot_t(int N) : x(N*N), y(N*N), n(N){
    generate(x.begin(), x.end(), gen_rand);
    generate(y.begin(), y.end(), gen_rand);
  }
};

void run_ddot(ddot_t d){
  cblas_ddot(d.n, &d.x.front(), 1, &d.y.front(), 1);
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

