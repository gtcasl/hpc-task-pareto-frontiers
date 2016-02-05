#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <lwperf/lwperf.h>
#include <omp.h>

#include <sys/time.h>
#include <sys/resource.h>
double mytimer(void) {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return tp.tv_sec + (double)tp.tv_usec/1e6;
}

extern void
generate_problem_27pt(
  int nx, int ny, int nz,
  int* nrows_ptr,
  int** nnzPerRow_ptr,
  int*** nonzerosInRow_ptr,
  double*** matrixRows_ptr,
  int** indices_ptr,
  double** matrix_ptr);

extern void
multiply(
  int nrows,
  double* vector,
  double* residual,
  double** matrixRows,
  int* nnzPerRow,
  int** nonzerosInRow);

int main(int argc, char** argv)
{
  if (argc < 4){
    std::cerr << "need 3 integers for nx,ny,nz on command line" << std::endl;
    return 1;
  }
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int nz = atoi(argv[3]);

  lwperf_t perf = lwperf_init("phi", "test", "database");
  lwperf_init_papi(perf);
  int nthreads = omp_get_max_threads();
  std::cout << "nthreads: " << nthreads << std::endl;
  lwperf_add_invariant(perf, "nthreads", nthreads);

  int nrows;
  double* matrix;
  int* indices;
  int** nonzerosInRow;
  double** matrixRows;
  int* nnzPerRow;

  generate_problem_27pt(nx,ny,nz,&nrows,&nnzPerRow,&nonzerosInRow,&matrixRows,&indices,&matrix);
  double* vector = new double[nrows];
  for (int r=0; r < nrows; ++r){
    vector[r] = 1.0;
  }

  double* residual = new double[nrows];

  double t_start = mytimer();
  int nloops = 50;
  lwperf_log(perf, "spmv");
  for (int l=0; l < nloops; l++){
    memset(residual,0,nrows*sizeof(double));
    multiply(nrows, vector, residual, matrixRows, nnzPerRow, nonzerosInRow);
  }
  lwperf_stop(perf, "spmv");
  double t_stop = mytimer();

  lwperf_finalize(perf);
  double t_total = t_stop - t_start;
  printf("For n=(%d,%d,%d) ran for %12.8fseconds\n", nx, ny, nz, t_total);

  delete[] matrix;
  delete[] indices;
  delete[] nonzerosInRow;
  delete[] matrixRows;
  delete[] nnzPerRow;

  delete[] vector;
  delete[] residual;
}


