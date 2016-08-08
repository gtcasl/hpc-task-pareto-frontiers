#ifdef no_mkl
#include <fake_mkl.h>
#else
#include <mkl.h>
#endif
#include <test.h>
#include <matrix.h>
#include "lu.h"

#if CHOLESKY_DEBUG
#define task_debug(...) printf(__VA_ARGS__)
#else
#define task_debug(...) 
#endif

extern __declspec(target(mic)) void copy(const double* src, double* dst, int size, int lda);

void
trsm_leftprofile(int k, int m, int size, int lda, double* L, double* A)
{
  //A overwritten with new U
  task_debug("Solving TRSM, L triangular  A(%d,%d) = L(%d,%d)*U(%d,%d)\n", 
             m, k, m, k, k, k);
#pragma offload target(mic:0) \
  in(L:length(0) alloc_if(0) free_if(0)) \
  in(A:length(0) alloc_if(0) free_if(0))
{
  double* Acopy = (double*) mkl_malloc(size * lda * sizeof(double), 64);
  copy(A, Acopy, size, lda);
  //solve AX = B
  //B overwrriten with X
  //char side = 'L';
  //char uplo = 'U';
  //char trans = 'N';
  //char diag = 'N';
  double alpha = 1.0;
  cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
    size, size, alpha,
    L, lda,
    Acopy, lda);
  mkl_free(Acopy);
}
}

void
trsm_rightprofile(int k, int m, int size, int lda, double* A, double* U)
{
  //A overwritten with new L 
  task_debug("Solving TRSM, U triangular  A(%d,%d) = L(%d,%d)*U(%d,%d)\n", 
             m, k, m, k, k, k);
#pragma offload target(mic:0) \
  in(U:length(0) alloc_if(0) free_if(0)) \
  in(A:length(0) alloc_if(0) free_if(0))
{
  double* Acopy = (double*) mkl_malloc(size * lda * sizeof(double), 64);
  copy(A, Acopy, size, lda);
  //solve AX = B
  //B overwrriten with X
  //char side = 'R';
  //char uplo = 'U';
  //char trans = 'N';
  //char diag = 'N';
  double alpha = 1.0;
  cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
    size, size, alpha,
    U, lda,
    Acopy, lda);
  mkl_free(Acopy);
}
}

void
getrfprofile(int k, int size, int lda, int* ipiv, double* A)
{
  task_debug("Running GETRF A(%d,%d)\n", k, k);
  //std::cout << A << std::endl;
  //for(int i = 0; i < size; i++){
  //  std::cout << A[i] << std::endl;
  //}
  int info;
#pragma offload target(mic:0) \
  in(A : length(0) alloc_if(0) free_if(0)) \
  in(ipiv : length(0) alloc_if(0) free_if(0)) \
  out(info)
{
  double* Acopy = (double*) mkl_malloc(size * lda * sizeof(double), 64);
  copy(A, Acopy, size, lda);
  info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, size, size, Acopy, lda, ipiv);
  mkl_free(Acopy);
}

  if (info != 0){
    fprintf(stderr, "FAILURE on DGETRF: %d\n", info);
    abort();
  }

}

void
lu_gemmprofile(
  int m, int n, int k,
  int size,
  int lda,
  bool ltrans, 
  bool rtrans,
  double alpha,
  double beta,
  double* product, 
  double* left, 
  double* right){
  task_debug("Adding  GEMM  L(%d,%d) += L(%d,%d)*L(%d,%d)\n", 
    m, n, m, k, k, n);
  CBLAS_TRANSPOSE lt = ltrans ? CblasNoTrans : CblasTrans; //GD fortran
  CBLAS_TRANSPOSE rt = rtrans ? CblasNoTrans : CblasTrans; //GD fortran
#pragma offload target(mic:0) \
  in(left : length(0) alloc_if(0) free_if(0)) \
  in(right : length(0) alloc_if(0) free_if(0)) \
  in(product : length(0) alloc_if(0) free_if(0))
{
  double* Pcopy = (double*) mkl_malloc(size * lda * sizeof(double), 64);
  copy(product, Pcopy, size, lda);
  cblas_dgemm(CblasColMajor, lt, rt, size, size, size,
    alpha, left, lda, 
    right, lda, 
    beta, Pcopy, lda);
  mkl_free(Pcopy);
}
}

class TaskMap {
 public:
  TaskMap(int size) : size_(size), tasks_(size*size, NULL){}
  
  Task*& operator()(int i, int j){
    int offset = i*size_ + j;
    return tasks_[offset];
  }
 private:
  std::vector<Task*> tasks_;
  int size_;
};

Task*
initProfilingDag(Matrix& A, int* ipiv_all)
{
  Task* root = 0;

  int nBlocks = A.blockGridSize();
  int blockSize = A.blockSize();
  int lda = nBlocks * blockSize;

  TaskMap gets(nBlocks);
  TaskMap trs(nBlocks);
  TaskMap gemms(nBlocks);
   
  for (int k=0; k < nBlocks; ++k){
    double* Akk = A.block(k,k);
    int* ipiv = ipiv_all + k*blockSize;
    Task* diagprofile = new_task(getrfprofile, k, blockSize, lda, ipiv, Akk);
    Task* diag = new_mutating_task(getrf, k, blockSize, lda, ipiv, Akk);
    diag->dependsOn(diagprofile);
    dep_debug("GETRF (%d,%d) = %p\n", k,k,diag);
    gets(k,k) = diag;
    if (k == 0){
      root = diagprofile;
    } else {
      //the potrf will only directly depend on the most recent syrk
      diagprofile->dependsOn(gemms(k,k));
    }
    for (int m=k+1; m < nBlocks; ++m){
      double* Amk = A.block(m,k);
      Task* solveprofile = new_task(trsm_leftprofile, k, m, blockSize, lda, Akk, Amk);
      Task* solve = new_mutating_task(trsm_left, k, m, blockSize, lda, Akk, Amk);
      dep_debug(" TRSM (%d,%d) = %p\n", m, k, solve);
      solve->dependsOn(solveprofile);
      trs(m,k) = solve;
      solveprofile->dependsOn(gets(k,k));
      if (m >= 2){
        solveprofile->dependsOn(gemms(m,k));
      }

      double* Akm = A.block(k,m);
      solveprofile = new_task(trsm_rightprofile, k, m, blockSize, lda, Akk, Akm);
      solve = new_mutating_task(trsm_right, k, m, blockSize, lda, Akk, Akm);
      dep_debug(" TRSM (%d,%d) = %p\n", k, m, solve);
      solve->dependsOn(solveprofile);
      trs(k,m) = solve;
      solveprofile->dependsOn(gets(k,k));
      if (m >= 2){
        solveprofile->dependsOn(gemms(k,m));
      }
    }
    for (int m=k+1; m < nBlocks; ++m){
      double* Amk = A.block(m,k);
      for (int n=k+1; n < nBlocks; ++n){
        double* Akn = A.block(k,n);
        double* Amn = A.block(m,n);
        Task* multprofile = new_task(lu_gemmprofile, m, n, k, blockSize, lda, true, false, -1.0, 1.0, Amn, Amk, Akn);
        Task* mult = new_mutating_task(lu_gemm, m, n, k, blockSize, lda, true, false, -1.0, 1.0, Amn, Amk, Akn);
        dep_debug(" GEMM (%d,%d) = %p\n", m,n,mult);
        mult->dependsOn(multprofile);
        gemms(m,n) = mult;
        multprofile->dependsOn(trs(m,k));
        multprofile->dependsOn(trs(k,n));
      }
    }
  }
  return root;
}

int luprofile(Scheduler* sch, int argc, char** argv)
{
  if(argc != 3){
    std::cerr << "Usage: " << argv[0] << " <nblocks> <blocksize>" << std::endl;
    return -1;
  }

  int nBlocks = atoi(argv[1]);
  int blockSize = atoi(argv[2]);
  auto nrows = nBlocks * blockSize;
  Matrix A(nBlocks, blockSize);
  Matrix L(nBlocks, blockSize);
  int* ipiv = (int*) malloc(nrows);

  std::cout << "Filling" << std::endl;
  A.randomFill();

#pragma offload_transfer target(mic:0) \
  in(A.storage_ : length(nrows * nrows) align(64) free_if(0) alloc_if(1)) \
  in(ipiv : length(nrows) align(64) free_if(0) alloc_if(1))

  Task* root = 0;
  root = initProfilingDag(A, ipiv);
  std::cout << "Running LU\n";

  sch->run(root);

  std::cout << "LU completed\n";
  fflush(stdout);

  //A.print("factorized A");


  return 0;
}

RegisterTest("luprofiling", luprofile);
