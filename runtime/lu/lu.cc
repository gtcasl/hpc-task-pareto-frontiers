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

void
trsm_left(int k, int m, int size, int lda, double* L, double* A)
{
  //A overwritten with new U
  task_debug("Solving TRSM, L triangular  A(%d,%d) = L(%d,%d)*U(%d,%d)\n", 
             m, k, m, k, k, k);
#pragma offload target(mic:0) \
  in(L:length(0) alloc_if(0) free_if(0)) \
  in(A:length(0) alloc_if(0) free_if(0))
{
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
    A, lda);
}
}

void
trsm_right(int k, int m, int size, int lda, double* A, double* U)
{
  //A overwritten with new L 
  task_debug("Solving TRSM, U triangular  A(%d,%d) = L(%d,%d)*U(%d,%d)\n", 
             m, k, m, k, k, k);
#pragma offload target(mic:0) \
  in(U:length(0) alloc_if(0) free_if(0)) \
  in(A:length(0) alloc_if(0) free_if(0))
{
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
    A, lda);
}
}

void
getrf(int k, int size, int lda, int* ipiv, double* A)
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
  info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, size, size, A, lda, ipiv);
}

/*
 * TODO: this is dumb. if this fails, it probably fails faster than a normal 
 * execution time, without touching everything it needs. It doesn't make much
 * practical sense, other than some handwaving to make it like we do the correct
 * number of flops.
  if (info != 0){
    fprintf(stderr, "FAILURE on DGETRF: %d\n", info);
    abort();
  }
  */

}

void
lu_gemm(
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
  task_debug("GEMM  L(%d,%d) += L(%d,%d)*L(%d,%d)\n", 
    m, n, m, k, k, n);
  CBLAS_TRANSPOSE lt = ltrans ? CblasNoTrans : CblasTrans; //GD fortran
  CBLAS_TRANSPOSE rt = rtrans ? CblasNoTrans : CblasTrans; //GD fortran
#pragma offload target(mic:0) \
  in(left : length(0) alloc_if(0) free_if(0)) \
  in(right : length(0) alloc_if(0) free_if(0)) \
  in(product : length(0) alloc_if(0) free_if(0))
{
  cblas_dgemm(CblasColMajor, lt, rt, size, size, size,
    alpha, left, lda, 
    right, lda, 
    beta, product, lda);
}
  task_debug("GEMM  L(%d,%d) += L(%d,%d)*L(%d,%d)\n", 
    m, n, m, k, k, n);
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
initDag(Matrix& A, int* ipiv_all)
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
    Task* diag = new_task(getrf, k, blockSize, lda, ipiv, Akk);
    dep_debug("GETRF (%d,%d) = %p\n", k,k,diag);
    gets(k,k) = diag;
    if (k == 0){
      root = diag;
    } else {
      //the potrf will only directly depend on the most recent syrk
      diag->dependsOn(gemms(k,k));
    }
    for (int m=k+1; m < nBlocks; ++m){
      double* Amk = A.block(m,k);
      Task* solve = new_task(trsm_left, k, m, blockSize, lda, Akk, Amk);
      dep_debug(" TRSM (%d,%d) = %p\n", m, k, solve);
      trs(m,k) = solve;
      solve->dependsOn(gets(k,k));
      if (m >= 2){
        solve->dependsOn(gemms(m,k));
      }

      double* Akm = A.block(k,m);
      solve = new_task(trsm_right, k, m, blockSize, lda, Akk, Akm);
      dep_debug(" TRSM (%d,%d) = %p\n", k, m, solve);
      trs(k,m) = solve;
      solve->dependsOn(gets(k,k));
      if (m >= 2){
        solve->dependsOn(gemms(k,m));
      }
    }
    for (int m=k+1; m < nBlocks; ++m){
      double* Amk = A.block(m,k);
      for (int n=k+1; n < nBlocks; ++n){
        double* Akn = A.block(k,n);
        double* Amn = A.block(m,n);
        Task* mult = new_task(lu_gemm, m, n, k, blockSize, lda, true, false, -1.0, 1.0, Amn, Amk, Akn);
        dep_debug(" GEMM (%d,%d) = %p\n", m,n,mult);
        gemms(m,n) = mult;
        mult->dependsOn(trs(m,k));
        mult->dependsOn(trs(k,n));
      }
    }
  }
  return root;
}

int lu(Scheduler* sch, int argc, char** argv)
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
  root = initDag(A, ipiv);
  std::cout << "Running LU\n";

  sch->run(root);

  std::cout << "LU completed\n";
  fflush(stdout);

  //A.print("factorized A");


  return 0;
}

RegisterTest("lu", lu);
