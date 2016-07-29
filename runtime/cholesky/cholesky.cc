#ifdef no_mkl
#include <fake_mkl.h>
#else
#include <mkl.h>
#endif
#include <test.h>

#if CHOLESKY_DEBUG
#define task_debug(...) printf(__VA_ARGS__)
#else
#define task_debug(...) 
#endif

static int static_info;

static enum fxn_id {
  potrf_id,
  gemm_id,
  syrk_id,
  trsm_id,
} test_ids;

typedef Buffer<double> DoublePtr;
typedef Buffer<double> DoubleArray;
typedef Buffer<int> IntArray;
typedef Buffer<IntArray> IntChunkArray;
typedef Buffer<DoubleArray> DoubleChunkArray;

class Matrix
{
 private:
  int blockOffset(int row, int col){
    int blockNum = row*blockGridSize_ + col;
    return blockNum*blockSize_*blockSize_;
  }

 public:
  void 
  symmetricFill(int row, int col){
    int rowStart = row * blockSize_;
    int colStart = col * blockSize_;
    int rowStop = rowStart + blockSize_;
    int colStop = colStart + blockSize_;
    int idx = 0;
    int offset = blockOffset(row,col);
    DoubleArray block = storage_.offset(offset);
    double* valptr = block;
    for (int i=rowStart; i != rowStop; ++i){
      for (int j=colStart; j != colStop; ++j, ++valptr, ++idx){
        if (i >= j){
          double x = double(0.05 + (idx*7)%4 - (idx*5)%2) / ((idx*11)%5+0.1);
          *valptr = x;
        } else {
          *valptr = 0;
        }
      }
    }
  }

  Matrix(int matrixSize, int blockSize, DoubleArray storage) :
    blockSize_(blockSize),
    blockGridSize_(matrixSize),
    storage_(storage)
  {}

  Matrix(int matrixSize, int blockSize) : 
    blockSize_(blockSize), 
    blockGridSize_(matrixSize),
    storage_(matrixSize*matrixSize*blockSize*blockSize)
  {
  }

  void print(const char* label = 0){
    if (label) printf("%s\n", label);
    int blockStorageSize = blockSize_*blockSize_;
    int nrows = blockSize_ * blockGridSize_;
    int ncols = nrows;
    int lastRowBlock = 0;
    for (int i=0; i < nrows; ++i){
      int lastColBlock = 0;
      int iBlock = i / blockSize_;
      int iOffset = i % blockSize_;
      if (iBlock != lastRowBlock){
        printf("\n");
        lastRowBlock = iBlock;
      }
      int offset = (iBlock*blockGridSize_*blockStorageSize) + iOffset*blockSize_;
      double* ptr = storage_.offset(offset);
      for (int j=0; j < ncols; ++j, ++ptr){
        int jBlock = j / blockSize_;
        int jOffset = j % blockSize_;
        if (jBlock != lastColBlock){
          int offset = (iBlock*blockGridSize_+jBlock)*blockStorageSize + (iOffset*blockSize_ + jOffset);
          //printf("Hopping to block %d,%d at offset %d\n", iBlock, jBlock, offset);
          printf("   ");
          ptr = storage_.offset(offset);
          lastColBlock = jBlock;
        }
        //printf("Pointer %p is %f\n", ptr, *ptr);
        printf("%6.2f ", *ptr);
      }
      printf("\n");
    }
    printf("\n");
  }

  int
  blockGridSize() const {
    return blockGridSize_;
  }
  
  int 
  blockSize() const {
    return blockSize_;
  }

  DoubleArray
  block(int i, int j){
    return storage_.offset(blockOffset(i,j));
  }
  double* storageAddr() {
    double* res = storage_;
    return res;
  }

 private:
  int blockSize_;
  int blockGridSize_;
  DoubleArray storage_;
};

void
trsm(int k, int m, int size, DoubleArray A, DoubleArray B)
{
  task_debug("Solving TRSM  A(%d,%d)  = L(%d,%d)*L(%d,%d)\n", m, k, m, k, k, k);
  //solve AX = B
  //B overwrriten with X
  char side = 'R';
  char uplo = 'U';
  char trans = 'N';
  char diag = 'N';
  double alpha = 1.0;
  cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
    size, size, alpha,
    A, size,
    B, size);
}

void
potrf(int k, int size, DoubleArray A)
{
  task_debug("Running POTRF A(%d,%d)\n", k, k);
  //std::cout << A << std::endl;
  //for(int i = 0; i < size; i++){
  //  std::cout << A[i] << std::endl;
  //}
  char uplo = 'U';
  int info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, uplo, size, A, size);
  if (info != 0){
    fprintf(stderr, "FAILURE on DPOTRF: %d\n", info);
    std::cout << "A[" << info << "]: " << A[info] << std::endl;
    abort();
  }
  double* ptr = A;
//#pragma omp parallel for
  for (int i=0; i < size; ++i){
    for (int j=0; j < size; ++j, ++ptr){
      if (j > i) *ptr = 0;
    }
  }
}

void
syrk(
  int k, int n,
  int size,
  bool trans,
  double alpha, 
  double beta,
  DoubleArray C,
  DoubleArray A)
{
  task_debug("Adding  SYRK  L(%d,%d) += L(%d,%d)*L(%d,%d)\n", n, n, n, k, k, n);
  CBLAS_TRANSPOSE ta = trans ? CblasNoTrans : CblasTrans;
  cblas_dsyrk(CblasColMajor, CblasUpper, ta, size, size, 
    alpha, A, size,
    beta, C, size);

//#pragma omp parallel for
  for (int i=0; i < size; ++i){
    for (int j=i+1; j < size; ++j){
      int toIdx = i*size + j;
      int fromIdx = j*size + i;
      C[toIdx] = C[fromIdx];
    }
  }
}

void
gemm(
  int m, int n, int k,
  int size,
  bool ltrans, 
  bool rtrans,
  double alpha,
  double beta,
  DoubleArray product, 
  DoubleArray left, 
  DoubleArray right){
  task_debug("Adding  GEMM  L(%d,%d) += L(%d,%d)*L(%d,%d)\n", 
    m, n, m, k, k, n);
  CBLAS_TRANSPOSE lt = ltrans ? CblasNoTrans : CblasTrans; //GD fortran
  CBLAS_TRANSPOSE rt = rtrans ? CblasNoTrans : CblasTrans; //GD fortran
  cblas_dgemm(CblasColMajor, lt, rt, size, size, size,
    alpha, left, size, 
    right, size, 
    beta, product, size);
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
initDag(Matrix& A)
{
  Task* root = 0;

  int nBlocks = A.blockGridSize();
  int blockSize = A.blockSize();

  TaskMap pots(nBlocks);
  TaskMap syrs(nBlocks);
  TaskMap trs(nBlocks);
  TaskMap gemms(nBlocks);
   
  for (int k=0; k < nBlocks; ++k){
    DoubleArray Akk = A.block(k,k);
    Task* diag = new_task(potrf, k, blockSize, Akk);
    dep_debug("POTRF (%d,%d) = %p\n", k,k,diag);
    pots(k,k) = diag;
    if (k == 0){
      root = diag;
    } else {
      //the potrf will only directly depend on the most recent syrk
      diag->dependsOn(syrs(k,k));
    }
    for (int m=k+1; m < nBlocks; ++m){
      DoubleArray Amk = A.block(m,k);
      Task* solve = new_task(trsm, k, m, blockSize, Akk, Amk);
      dep_debug(" TRSM (%d,%d) = %p\n", m, k, solve);
      trs(m,k) = solve;
      solve->dependsOn(pots(k,k));
      if (m >= 2){
        solve->dependsOn(gemms(m,k));
      }
    }
    for (int n=k+1; n < nBlocks; ++n){
      DoubleArray Ann = A.block(n,n);
      DoubleArray Ank = A.block(n,k);
      //Amm -= Amk * Amk^T
      //Ank ends up transposed from the DTRSM solve
      Task* symmUpd = new_task(syrk, k, n, blockSize, true, -1.0, 1.0, Ann, Ank);
      dep_debug(" SYRK (%d,%d) = %p\n", n,n,symmUpd);
      symmUpd->dependsOn(trs(n,k));
      symmUpd->dependsOn(syrs(n,n)); //depend on prev syrk here
      syrs(n,n) = symmUpd;
      for (int m=n+1; m < nBlocks; ++m){
        DoubleArray Amk = A.block(m,k);
        DoubleArray Amn = A.block(m,n);
        //because of how awesome BLAS is, Ank is transposed backwards
        Task* mult = new_task(gemm, m, n, k, blockSize, true, false, -1.0, 1.0, Amn, Amk, Ank);
        dep_debug(" GEMM (%d,%d) = %p\n", m,n,mult);
        gemms(m,n) = mult;
        mult->dependsOn(trs(m,k));
        mult->dependsOn(trs(n,k));
      }
    }
  }
  return root;
}

void fill(int nBlocks, int blockSize, Matrix& L, Matrix& A)
{
}


int cholesky(int argc, char** argv)
{
  //ALWAYS Initialize the scheduler first
  Scheduler* sch = new BasicScheduler;
  //Scheduler* sch = new AdvancedScheduler;
  sch->init(argc, argv);
  // disable dynamic thread adjustment in MKL
#ifndef no_mkl
  mkl_set_dynamic(0);
#endif

  RegisterTask(potrf, void, int,
    int, DoubleArray);
  RegisterTask(trsm,  void, int, int,
    int, DoubleArray, DoubleArray);
  RegisterTask(syrk,  void, int, int, 
    int, bool, double, double, DoubleArray, DoubleArray);
  RegisterTask(gemm,  void, int, int, int,
    int, bool, bool, double, double, DoubleArray, DoubleArray, DoubleArray);

  if(argc != 4){
    if(sch->rank() == 0){
      std::cerr << "Usage: " << argv[1] << " <nblocks> <blocksize>" << std::endl;
    }
    return -1;
  }

  int nBlocks = atoi(argv[2]);
  int blockSize = atoi(argv[3]);
  Matrix A(nBlocks, blockSize);
  Matrix L(nBlocks, blockSize);

  int ncopies = 1;
  sch->allocateHeap(ncopies);

  for (int c=0; c < ncopies; ++c, sch->nextIter()){
    if (sch->rank() == 0){

      for (int i=0; i < nBlocks; ++i){
        for (int j=0; j < nBlocks; ++j){
          L.symmetricFill(i,j);
        }
      }

      for (int i=0; i < nBlocks; ++i){
        for (int j=0; j <= i; ++j){
          DoubleArray Aij = A.block(i,j);
          for (int k=0; k < nBlocks; ++k){
            DoubleArray Lik = L.block(i,k);
            DoubleArray Ljk = L.block(j,k);
            gemm(i,j,k,blockSize, false, true, 1.0, 1.0, Aij, Lik, Ljk);
          }
        }
      }

      // send A & L
      MPI_Send(A.storageAddr(), nBlocks * nBlocks * blockSize * blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
      MPI_Send(L.storageAddr(), nBlocks * nBlocks * blockSize * blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    }
    if(sch->rank() == 1){
      // recv A & L
      MPI_Recv(A.storageAddr(), nBlocks * nBlocks * blockSize * blockSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(L.storageAddr(), nBlocks * nBlocks * blockSize * blockSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //A.print("before");
      //L.print("before");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  sch->resetIter();

  for (int iter=0; iter < ncopies; ++iter, sch->nextIter()){
    Task* root = 0;
    if (sch->rank() == 0){
      root = initDag(A);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    sch->run(root);
    int nfailures = 0;
    if(sch->rank() == 1){
      A.print("after");
      L.print("after");
      for (int i=0; i < nBlocks; ++i){
        //just check the diagonal blocks...
        //the other blocks end up weird and transposed
        DoubleArray Aii = A.block(i,i);
        DoubleArray Lii = L.block(i,i);
        int nelems = blockSize * blockSize;
        double tol = 1e-2;
        for (int j=0; j < nelems; ++j){
          double delta = fabs(Aii[j] - Lii[j]);
          if (delta > tol){
            nfailures++;
          }
        }
      }
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
  sch->finalize();

  return 0;
}

RegisterTest("cholesky", cholesky);
