#ifdef no_mkl
#include <fake_mkl.h>
#else
#include <mkl.h>
#endif
#include <test.h>

#include "cholesky.h"

//#define CHOLESKY_DEBUG 1
#if CHOLESKY_DEBUG
#define task_debug(...) printf(__VA_ARGS__)
#else
#define task_debug(...) 
#endif

void
trsm(int k, int m, int size, int lda, double* A, double* B)
{
  task_debug("TRSM  A(%d,%d)  = L(%d,%d)*L(%d,%d)\n", m, k, m, k, k, k);
#pragma offload target(mic:0) \
  in(size) \
  in(lda) \
  in(A:length(0) alloc_if(0) free_if(0)) \
  in(B:length(0) alloc_if(0) free_if(0))
  {
  // A = Akk
  // B = Amk
  //solve AX = B
  //B overwrriten with X
  //char side = 'R';
  //char uplo = 'U';
  //char trans = 'N';
  //char diag = 'N';
  double alpha = 1.0;
  cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
    size, size, alpha,
    A, lda,
    B, lda);
  }
  task_debug("TRSM  A(%d,%d)  = L(%d,%d)*L(%d,%d)\n", m, k, m, k, k, k);
}

void
potrf(int k, int size, int lda, double* A)
{
  task_debug("POTRF A(%d,%d)\n", k, k);
  char uplo = 'L';
  int info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, uplo, size, A, lda);
  if (info != 0){
    fprintf(stderr, "FAILURE on DPOTRF: %d\n", info);
    std::cout << "A[" << info << "]: " << A[info] << std::endl;
    abort();
  }
  task_debug("POTRF A(%d,%d)\n", k, k);
}

void
syrk(
  int k, int n,
  int size,
  int lda,
  double* C,
  double* A)
{
  // C=Ann
  // A=Ank
  // C = C - A*A'
  task_debug("SYRK  L(%d,%d) += L(%d,%d)*L(%d,%d)\n", n, n, n, k, n, n);
  auto alpha = -1.0;
  auto beta = 1.0;
  cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, size, size, 
    alpha, A, lda,
    beta, C, lda);
  task_debug("SYRK  L(%d,%d) += L(%d,%d)*L(%d,%d)\n", n, n, n, k, n, n);
}

void
gemm(
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
  cblas_dgemm(CblasColMajor, lt, rt, size, size, size,
    alpha, left, lda, 
    right, lda, 
    beta, product, lda);
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
  int size_;
  std::vector<Task*> tasks_;
};

Task*
initDag(Matrix& A)
{
  Task* root = 0;

  int nBlocks = A.blockGridSize();
  int blockSize = A.blockSize();
  int lda = nBlocks * blockSize;

  TaskMap pots(nBlocks);
  TaskMap syrs(nBlocks);
  TaskMap trs(nBlocks);
  TaskMap gemms(nBlocks);
   
  for (int k=0; k < nBlocks; ++k){
    double* Akk = A.block(k,k);
    Task* diag = new_task(potrf, k, blockSize, lda, Akk);
    dep_debug("POTRF (%d,%d) = %p\n", k,k,diag);
    pots(k,k) = diag;
    if (k == 0){
      root = diag;
    } else {
      //the potrf will only directly depend on the most recent syrk
      diag->dependsOn(syrs(k,k));
    }
    for (int m=k+1; m < nBlocks; ++m){
      double* Amk = A.block(m,k);
      Task* solve = new_task(trsm, k, m, blockSize, lda, Akk, Amk);
      dep_debug(" TRSM (%d,%d) = %p\n", m, k, solve);
      trs(m,k) = solve;
      solve->dependsOn(pots(k,k));
      if (m >= 2){
        solve->dependsOn(gemms(m,k));
      }
    }
    for (int n=k+1; n < nBlocks; ++n){
      double* Ann = A.block(n,n);
      double* Ank = A.block(n,k);
      //Amm -= Amk * Amk^T
      //Ank ends up transposed from the DTRSM solve
      Task* symmUpd = new_task(syrk, k, n, blockSize, lda, Ann, Ank);
      dep_debug(" SYRK (%d,%d) = %p\n", n,n,symmUpd);
      symmUpd->dependsOn(trs(n,k));
      symmUpd->dependsOn(syrs(n,n)); //depend on prev syrk here
      syrs(n,n) = symmUpd;
      for (int m=n+1; m < nBlocks; ++m){
        double* Amk = A.block(m,k);
        double* Amn = A.block(m,n);
        //because of how awesome BLAS is, Ank is transposed backwards
        Task* mult = new_task(gemm, m, n, k, blockSize, lda, true, false, -1.0, 1.0, Amn, Amk, Ank);
        dep_debug(" GEMM (%d,%d) = %p\n", m,n,mult);
        gemms(m,n) = mult;
        mult->dependsOn(trs(m,k));
        mult->dependsOn(trs(n,k));
      }
    }
  }
  return root;
}

int cholesky(Scheduler* sch, int argc, char** argv)
{
  // disable dynamic thread adjustment in MKL
#ifndef no_mkl
  mkl_set_dynamic(0);
#endif

  RegisterTask(potrf, void, int, int,
    int, double*);
  RegisterTask(trsm,  void, int, int, int,
    int, double*, double*);
  RegisterTask(syrk,  void, int, int, int,
    int, double*, double*);
  RegisterTask(gemm,  void, int, int, int, int,
    int, bool, bool, double, double, double*, double*, double*);

  if(argc != 3){
    std::cerr << "Usage: " << argv[0] << " <nblocks> <blocksize>" << std::endl;
    return -1;
  }

  std::cout << "Initializing cholesky" << std::endl;

  int nBlocks = atoi(argv[1]);
  int blockSize = atoi(argv[2]);
  auto nrows = nBlocks * blockSize;
  Matrix A(nBlocks, blockSize);
  Matrix L(nBlocks, blockSize);

  int ncopies = 1;

  for (int c=0; c < ncopies; ++c, sch->nextIter()){
#ifdef _OPENMP
    omp_set_num_threads(20);
#endif
    std::cout << "Filling" << std::endl;
    L.symmetricFill();
    A.randomFill();
    std::cout << "Blas" << std::endl;
    // tmp = A' x L
    double* tmp = (double*) malloc(nrows * nrows * sizeof(double));
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nrows, nrows, nrows,
                1.0, A.storageAddr(), nrows, L.storageAddr(), nrows, 0.0,
                tmp, nrows);
    // L = tmp x A
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nrows, nrows, nrows,
                1.0, tmp, nrows, A.storageAddr(), nrows, 0.0,
                L.storageAddr(), nrows);
    free(tmp);
    // A = copy(L)
    cblas_dcopy(nrows * nrows, L.storageAddr(), 1, A.storageAddr(), 1);
    std::cout << "Blas done" << std::endl;
  }
  // return 0;

#pragma offload_transfer target(mic:0) \
  in(L.storage_ : length(nrows * nrows) align(64) free_if(0) alloc_if(1))

  for (int iter=0; iter < ncopies; ++iter, sch->nextIter()){
    Task* root = 0;
    root = initDag(L);
    std::cout << "Running Cholesky\n";

    sch->run(root);
    int nfailures = 0;
    std::cout << "Performing validation" << std::endl;
    /*
#ifdef _OPENMP
    omp_set_num_threads(20);
#endif
    for (int i=0; i < nrows; ++i){
      for (int j=i+1; j < nrows; ++j){
        int toIdx = j*nrows + i;
        L.storageAddr()[toIdx] = 0.0;
      }
    }
    std::cout << "blas" << std::endl;
    double* tmp = (double*) malloc(nrows * nrows * sizeof(double));
    cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, nrows, nrows,
                1.0, L.storageAddr(), nrows, 0.0, tmp, nrows);
    std::cout << "test" << std::endl;
#pragma omp parallel for
    for(int i = 0; i < nrows; i++){
      for(int j = 0; j < i+1; j++){
        double tol = 1e-4;
        auto gold_val = A.storageAddr()[j*nrows + i];
        auto test_val = tmp[j*nrows + i];
        if(fabs(gold_val - test_val) > tol){
          nfailures++;
        }
      }
    }
    free(tmp);
    if (nfailures){
      printf("Cholesky failed with %d wrong elements on iteration %d\n",
        nfailures, iter);
    } else {
      printf("Cholesky passed validation test on iteration %d\n", iter);
    }
  }
  */
  }

  fflush(stdout);

  return 0;
}

RegisterTest("cholesky", cholesky);
