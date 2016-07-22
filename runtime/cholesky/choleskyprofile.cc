#ifdef no_mkl
#include <fake_mkl.h>
#else
#include <mkl.h>
#endif
#include <test.h>

#include "cholesky.h"

#if CHOLESKY_DEBUG
#define task_debug(...) printf(__VA_ARGS__)
#else
#define task_debug(...) 
#endif

void copy(const double* src, double* dst, int size, int lda){
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      dst[i*lda + j] = src[i*lda +j];
    }
  }
}

void dotrsm(int k, int m, int size, int lda, DoubleArray A, DoubleArray B){
  trsm(k, m, size, lda, A, B);
}
void dopotrf(int k, int size, int lda, DoubleArray A){
  potrf(k, size, lda, A);
}
void dosyrk(int k, int n, int size, int lda, DoubleArray C, DoubleArray A){
  syrk(k, n, size, lda, C, A);
}
void dogemm(int m, int n, int k, int size, int lda, bool ltrans,  bool rtrans, double alpha, double beta, DoubleArray product,  DoubleArray left,  DoubleArray right){
  gemm(m, n, k, size, lda, ltrans,  rtrans, alpha, beta, product,  left,  right);
}

void
trsmprofile(int k, int m, int size, int lda, DoubleArray A, DoubleArray B)
{
  task_debug("Solving TRSM  A(%d,%d)  = L(%d,%d)*L(%d,%d)\n", m, k, m, k, k, k);
  double* Bcopy = (double*) malloc(size * lda * sizeof(double));
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

  free(Bcopy);
}

void
potrfprofile(int k, int size, int lda, DoubleArray A)
{
  double* Acopy = (double*) malloc(size * lda * sizeof(double));
  copy(A, Acopy, size, lda);
  task_debug("Running POTRF A(%d,%d)\n", k, k);
  char uplo = 'L';
  int info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, uplo, size, Acopy, lda);
  if (info != 0){
    fprintf(stderr, "FAILURE on DPOTRF: %d\n", info);
    std::cout << "A[" << info << "]: " << Acopy[info] << std::endl;
    abort();
  }
  free(Acopy);
}

void
syrkprofile(
  int k, int n,
  int size,
  int lda,
  DoubleArray C,
  DoubleArray A)
{
  double* Ccopy = (double*) malloc(size * lda * sizeof(double));
  copy(C, Ccopy, size, lda);
  task_debug("Adding  SYRK  L(%d,%d) += L(%d,%d)*L(%d,%d)\n", n, n, n, k, k, n);
  auto alpha = -1.0;
  auto beta = 1.0;
  cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, size, size, 
    alpha, A, lda,
    beta, Ccopy, lda);

  free(Ccopy);
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
  DoubleArray product, 
  DoubleArray left, 
  DoubleArray right){
  double* Pcopy = (double*) malloc(size * lda * sizeof(double));
  copy(product, Pcopy, size, lda);
  task_debug("Adding  GEMM  L(%d,%d) += L(%d,%d)*L(%d,%d)\n", 
    m, n, m, k, k, n);
  CBLAS_TRANSPOSE lt = ltrans ? CblasNoTrans : CblasTrans; //GD fortran
  CBLAS_TRANSPOSE rt = rtrans ? CblasNoTrans : CblasTrans; //GD fortran
  cblas_dgemm(CblasColMajor, lt, rt, size, size, size,
    alpha, left, lda, 
    right, lda, 
    beta, Pcopy, lda);
  free(Pcopy);
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
    DoubleArray Akk = A.block(k,k);
    Task* diagprofile = new_task(potrfprofile, k, blockSize, lda, Akk);
    Task* diag = new_task(dopotrf, k, blockSize, lda, Akk);
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
      DoubleArray Amk = A.block(m,k);
      Task* solveprofile = new_task(trsmprofile, k, m, blockSize, lda, Akk, Amk);
      Task* solve = new_task(dotrsm, k, m, blockSize, lda, Akk, Amk);
      dep_debug(" TRSM (%d,%d) = %p\n", m, k, solve);
      solve->dependsOn(solveprofile);
      trs(m,k) = solve;
      solveprofile->dependsOn(pots(k,k));
      if (m >= 2){
        solveprofile->dependsOn(gemms(m,k));
      }
    }
    for (int n=k+1; n < nBlocks; ++n){
      DoubleArray Ann = A.block(n,n);
      DoubleArray Ank = A.block(n,k);
      //Amm -= Amk * Amk^T
      //Ank ends up transposed from the DTRSM solve
      Task* symmUpdprofile = new_task(syrkprofile, k, n, blockSize, lda, Ann, Ank);
      Task* symmUpd = new_task(dosyrk, k, n, blockSize, lda, Ann, Ank);
      symmUpd->dependsOn(symmUpdprofile);
      dep_debug(" SYRK (%d,%d) = %p\n", n,n,symmUpd);
      symmUpdprofile->dependsOn(trs(n,k));
      symmUpdprofile->dependsOn(syrs(n,n)); //depend on prev syrk here
      syrs(n,n) = symmUpd;
      for (int m=n+1; m < nBlocks; ++m){
        DoubleArray Amk = A.block(m,k);
        DoubleArray Amn = A.block(m,n);
        //because of how awesome BLAS is, Ank is transposed backwards
        Task* multprofile = new_task(gemmprofile, m, n, k, blockSize, lda, true, false, -1.0, 1.0, Amn, Amk, Ank);
        Task* mult = new_task(dogemm, m, n, k, blockSize, lda, true, false, -1.0, 1.0, Amn, Amk, Ank);
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
  // disable dynamic thread adjustment in MKL
#ifndef no_mkl
  mkl_set_dynamic(0);
#endif

  RegisterTask(potrfprofile, void, int, int,
    int, DoubleArray);
  RegisterTask(trsmprofile,  void, int, int, int,
    int, DoubleArray, DoubleArray);
  RegisterTask(syrkprofile,  void, int, int, int,
    int, DoubleArray, DoubleArray);
  RegisterTask(gemmprofile,  void, int, int, int, int,
    int, bool, bool, double, double, DoubleArray, DoubleArray, DoubleArray);
  RegisterMutatingTask(dopotrf, void, int, int,
    int, DoubleArray);
  RegisterMutatingTask(dotrsm,  void, int, int, int,
    int, DoubleArray, DoubleArray);
  RegisterMutatingTask(dosyrk,  void, int, int, int,
    int, DoubleArray, DoubleArray);
  RegisterMutatingTask(dogemm,  void, int, int, int, int,
    int, bool, bool, double, double, DoubleArray, DoubleArray, DoubleArray);

  if(argc != 3){
    if(sch->rank() == 0){
      std::cerr << "Usage: " << argv[0] << " <nblocks> <blocksize>" << std::endl;
    }
    return -1;
  }

  int nBlocks = atoi(argv[1]);
  int blockSize = atoi(argv[2]);
  auto nrows = nBlocks * blockSize;
  Matrix A(nBlocks, blockSize);
  Matrix L(nBlocks, blockSize);

  int ncopies = 1;
  sch->allocateHeap(ncopies);

  for (int c=0; c < ncopies; ++c, sch->nextIter()){
    if (sch->rank() == 0){

      L.symmetricFill();
      A.randomFill();
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

      // send A & L
      MPI_Send(A.storageAddr(), nBlocks * nBlocks * blockSize * blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
      MPI_Send(L.storageAddr(), nBlocks * nBlocks * blockSize * blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    }
    if(sch->rank() == 1){
      // recv A & L
      MPI_Recv(A.storageAddr(), nBlocks * nBlocks * blockSize * blockSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(L.storageAddr(), nBlocks * nBlocks * blockSize * blockSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  sch->resetIter();

  for (int iter=0; iter < ncopies; ++iter, sch->nextIter()){
    Task* root = 0;
    if (sch->rank() == 0){
      root = initProfilingDag(L);
      std::cout << "Running Cholesky profiling\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    sch->run(root);
    int nfailures = 0;
    if(sch->rank() == 1){
      std::cout << "Performing validation\n";
      for (int i=0; i < nrows; ++i){
        for (int j=i+1; j < nrows; ++j){
          int toIdx = j*nrows + i;
          L.storageAddr()[toIdx] = 0.0;
        }
      }
      double* tmp = (double*) malloc(nrows * nrows * sizeof(double));
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nrows, nrows, nrows,
                  1.0, L.storageAddr(), nrows, L.storageAddr(), nrows,
                  0.0, tmp, nrows);
      for(int i = 0; i < nrows * nrows; i++){
        double tol = 1e-4;
        auto gold_val = A.storageAddr()[i];
        auto test_val = tmp[i];
        if(fabs(gold_val - test_val) > tol){
          nfailures++;
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
  }

  fflush(stdout);

  //A.print("factorized A");


  sch->deallocateHeap();

  return 0;
}

RegisterTest("choleskyprofiling", choleskyprofiling);
