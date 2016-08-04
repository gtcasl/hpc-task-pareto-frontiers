#include <mkl.h>
#include <test.h>

#include "cholesky.h"

#if CHOLESKY_DEBUG
#define task_debug(...) printf(__VA_ARGS__)
#else
#define task_debug(...) 
#endif

__declspec(target(mic)) void copy(const double* src, double* dst, int size, int lda){
#pragma omp parallel for
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      dst[i*lda + j] = src[i*lda +j];
    }
  }
}

void dotrsm(int k, int m, int size, int lda, double* A, double* B){
  trsm(k, m, size, lda, A, B);
}
void dopotrf(int k, int size, int lda, double* A){
  potrf(k, size, lda, A);
}
void dosyrk(int k, int n, int size, int lda, double* C, double* A){
  syrk(k, n, size, lda, C, A);
}
void dogemm(int m, int n, int k, int size, int lda, bool ltrans,  bool rtrans, double alpha, double beta, double* product,  double* left,  double* right){
  gemm(m, n, k, size, lda, ltrans,  rtrans, alpha, beta, product,  left,  right);
}

void
trsmprofile(int k, int m, int size, int lda, double* A, double* B)
{
  task_debug("Solving TRSM  A(%d,%d)  = L(%d,%d)*L(%d,%d)\n", m, k, m, k, k, k);
#pragma offload target(mic:0) \
  in(A:length(0) alloc_if(0) free_if(0)) \
  in(B:length(0) alloc_if(0) free_if(0))
{
  double* Bcopy = (double*) mkl_malloc(size * lda * sizeof(double), 64);
  copy(B, Bcopy, size, lda);
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
    Bcopy, lda);

  mkl_free(Bcopy);
}
}

void
potrfprofile(int k, int size, int lda, double* A)
{
  task_debug("Running POTRF A(%d,%d)\n", k, k);
  char uplo = 'L';
  int info;
#pragma offload target(mic:0) \
  in(A : length(0) alloc_if(0) free_if(0)) \
  out(info)
{
  double* Acopy = (double*) mkl_malloc(size * lda * sizeof(double), 64);
  copy(A, Acopy, size, lda);
  info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, uplo, size, Acopy, lda);
  mkl_free(Acopy);
}
  if (info != 0){
    fprintf(stderr, "FAILURE on DPOTRF: %d\n", info);
    abort();
  }
}

void
syrkprofile(
  int k, int n,
  int size,
  int lda,
  double* C,
  double* A)
{
  task_debug("Adding  SYRK  L(%d,%d) += L(%d,%d)*L(%d,%d)\n", n, n, n, k, k, n);
  auto alpha = -1.0;
  auto beta = 1.0;
#pragma offload target(mic:0) \
  in(A : length(0) alloc_if(0) free_if(0)) \
  in(C : length(0) alloc_if(0) free_if(0)) 
{
  double* Ccopy = (double*) mkl_malloc(size * lda * sizeof(double), 64);
  copy(C, Ccopy, size, lda);
  cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, size, size, 
    alpha, A, lda,
    beta, Ccopy, lda);

  mkl_free(Ccopy);
}
}

void
gemmprofile(
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
  int size_;
  std::vector<Task*> tasks_;
};

Task*
initProfilingDag(Matrix& A)
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
    Task* diagprofile = new_task(potrfprofile, k, blockSize, lda, Akk);
    Task* diag = new_mutating_task(dopotrf, k, blockSize, lda, Akk);
    dep_debug("POTRF (%d,%d) = %p\n", k,k,diag);
    diag->dependsOn(diagprofile);
    pots(k,k) = diag;
    if (k == 0){
      root = diagprofile;
    } else {
      //the potrf will only directly depend on the most recent syrk
      diagprofile->dependsOn(syrs(k,k));
    }
    for (int m=k+1; m < nBlocks; ++m){
      double* Amk = A.block(m,k);
      Task* solveprofile = new_task(trsmprofile, k, m, blockSize, lda, Akk, Amk);
      Task* solve = new_mutating_task(dotrsm, k, m, blockSize, lda, Akk, Amk);
      dep_debug(" TRSM (%d,%d) = %p\n", m, k, solve);
      solve->dependsOn(solveprofile);
      trs(m,k) = solve;
      solveprofile->dependsOn(pots(k,k));
      if (m >= 2){
        solveprofile->dependsOn(gemms(m,k));
      }
    }
    for (int n=k+1; n < nBlocks; ++n){
      double* Ann = A.block(n,n);
      double* Ank = A.block(n,k);
      //Amm -= Amk * Amk^T
      //Ank ends up transposed from the DTRSM solve
      Task* symmUpdprofile = new_task(syrkprofile, k, n, blockSize, lda, Ann, Ank);
      Task* symmUpd = new_mutating_task(dosyrk, k, n, blockSize, lda, Ann, Ank);
      symmUpd->dependsOn(symmUpdprofile);
      dep_debug(" SYRK (%d,%d) = %p\n", n,n,symmUpd);
      symmUpdprofile->dependsOn(trs(n,k));
      symmUpdprofile->dependsOn(syrs(n,n)); //depend on prev syrk here
      syrs(n,n) = symmUpd;
      for (int m=n+1; m < nBlocks; ++m){
        double* Amk = A.block(m,k);
        double* Amn = A.block(m,n);
        //because of how awesome BLAS is, Ank is transposed backwards
        Task* multprofile = new_task(gemmprofile, m, n, k, blockSize, lda, true, false, -1.0, 1.0, Amn, Amk, Ank);
        Task* mult = new_mutating_task(dogemm, m, n, k, blockSize, lda, true, false, -1.0, 1.0, Amn, Amk, Ank);
        mult->dependsOn(multprofile);
        dep_debug(" GEMM (%d,%d) = %p\n", m,n,mult);
        gemms(m,n) = mult;
        multprofile->dependsOn(trs(m,k));
        multprofile->dependsOn(trs(n,k));
      }
    }
  }
  return root;
}

int choleskyprofiling(Scheduler* sch, int argc, char** argv)
{
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
  // A = copy(L)
  cblas_dcopy(nrows * nrows, L.storageAddr(), 1, A.storageAddr(), 1);
  std::cout << "Blas done" << std::endl;
  // return 0;

#pragma offload_transfer target(mic:0) \
  in(L.storage_ : length(nrows * nrows) align(64) free_if(0) alloc_if(1))

  Task* root = 0;
  root = initProfilingDag(L);
  std::cout << "Running Cholesky\n";

  sch->run(root);
  int nfailures = 0;
  std::cout << "Performing validation" << std::endl;
#pragma offload_transfer target(mic:0) \
  out(L.storage_ : length(nrows * nrows) free_if(1))

  for (int i=0; i < nrows; ++i){
    for (int j=i+1; j < nrows; ++j){
      int toIdx = j*nrows + i;
      L.storageAddr()[toIdx] = 0.0;
    }
  }
  std::cout << "blas" << std::endl;
  cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, nrows, nrows,
              1.0, L.storageAddr(), nrows, 0.0, tmp, nrows);
  std::cout << "test" << std::endl;
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
    printf("Cholesky failed with %d wrong elements\n",
      nfailures);
  } else {
    printf("Cholesky passed validation test\n");
  }

  fflush(stdout);

  return 0;
}

RegisterTest("choleskyprofiling", choleskyprofiling);
