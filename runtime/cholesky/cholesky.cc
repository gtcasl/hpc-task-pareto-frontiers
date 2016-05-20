#include <test.h>

static enum fxn_id {
  potrf_id,
  gemm_id,
  syrk_id,
  trsm_id,
} test_ids;

extern "C" 
int dgemm(char *transa, char *transb, int *m, int *n, int *k, 
  double *alpha, double *a, int *lda, 
  double *b, int *ldb, double *beta, 
  double *c, int *ldc);

extern "C"
int dsyrk(char *uplo, char *trans, int *n, int *k, 
  double *alpha, double *a, int *lda, double *beta, 
  double *c, int *ldc);

extern "C" 
int dpotrf(char *uplo, int *n, double *a, int* lda, int* info);

extern "C" 
int dtrsm(char *side, char *uplo, char *transa, char *diag, 
  int *m, int *n, 
  double* alpha, double* a, int* lda, 
  double* b, int* ldb);

typedef Buffer<double> DoublePtr;
typedef Buffer<double> DoubleArray;
typedef Buffer<int> IntArray;
typedef Buffer<IntArray> IntChunkArray;
typedef Buffer<DoubleArray> DoubleChunkArray;

void
trsm(int size, DoubleArray A, DoubleArray B)
{
  //solve AX = B
  //B overwrriten with X
  char side = 'R';
  char uplo = 'U';
  char trans = 'N';
  char diag = 'N';
  double alpha = 1.0;
  dtrsm(&side, &uplo, &trans, &diag,
    &size, &size, &alpha,
    A.buffer, &size,
    B.buffer, &size);
}

void
potrf(int size, DoubleArray A)
{
  char uplo = 'U';
  int info;
  dpotrf(&uplo, &size, A.buffer, &size, &info);
  if (info != 0){
    fprintf(stderr, "FAILURE on DPOTRF: %d\n", info);
    abort();
  }
  double* ptr = A.buffer;
#pragma omp parallel for
  for (int i=0; i < size; ++i){
    for (int j=0; j < size; ++j, ++ptr){
      if (j > i) *ptr = 0;
    }
  }
}

void
syrk(
  int size,
  bool trans,
  double alpha, 
  double beta,
  DoubleArray C,
  DoubleArray A)
{
  char ta = trans ? 'N' : 'T';
  char uplo = 'U';
  dsyrk(&uplo, &ta, &size, &size, 
    &alpha, A.buffer, &size,
    &beta, C.buffer, &size);

#pragma omp parallel for
  for (int i=0; i < size; ++i){
    for (int j=i+1; j < size; ++j){
      int toIdx = i*size + j;
      int fromIdx = j*size + i;
      C[toIdx] = C[fromIdx];
    }
  }
}

void
gemm(int size,
  bool ltrans, 
  bool rtrans,
  double alpha,
  double beta,
  DoubleArray product, 
  DoubleArray left, 
  DoubleArray right){
  char lt = ltrans ? 'N' : 'T'; //GD fortran
  char rt = rtrans ? 'N' : 'T'; //GD fortran
  dgemm(&lt, &rt, &size, &size, &size,
    &alpha, left.buffer, &size, 
    right.buffer, &size, 
    &beta, product.buffer, &size);
}

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
    DoubleArray block = storage_.offset(blockOffset(row,col));
    double* valptr = block.buffer;
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
      double* ptr = storage_.offset(offset).buffer;
      for (int j=0; j < ncols; ++j, ++ptr){
        int jBlock = j / blockSize_;
        int jOffset = j % blockSize_;
        if (jBlock != lastColBlock){
          int offset = (iBlock*blockGridSize_+jBlock)*blockStorageSize + (iOffset*blockSize_ + jOffset);
          printf("    ");
          ptr = storage_.offset(offset).buffer;
          lastColBlock = jBlock;
        }
        //printf("Pointer %p is %f\n", ptr, *ptr);
        printf("%10.6f", *ptr);
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

 private:
  int blockSize_;
  int blockGridSize_;
  DoubleArray storage_;
};

Task*
initDag(Matrix& A)
{
  Task* root = 0;

  int nBlocks = A.blockGridSize();
  int blockSize = A.blockSize();
   
  for (int k=0; k < nBlocks; ++k){
    DoubleArray Akk = A.block(k,k);
    Task* diag = new_task(potrf, blockSize, Akk);
    for (int m=k+1; m < nBlocks; ++m){
      DoubleArray Amk = A.block(m,k);
      Task* solve = new_task(trsm, blockSize, Akk, Amk);
    }
    for (int n=k+1; n < nBlocks; ++n){
      DoubleArray Ann = A.block(n,n);
      DoubleArray Ank = A.block(n,k);
      //Amm -= Amk * Amk^T
      //Ank ends up transposed from the DTRSM solve
      new_task(syrk, blockSize, true, -1.0, 1.0, Ann, Ank);
      for (int m=n+1; m < nBlocks; ++m){
        DoubleArray Amk = A.block(m,k);
        DoubleArray Amn = A.block(m,n);
        //because of how awesome BLAS is, Ank is transposed backwards
        new_task(gemm, blockSize, true, false, -1.0, 1.0, Amn, Amk, Ank);
      }
    }
  }
  return root;
}


int cholesky(int argc, char** argv)
{
  //ALWAYS Initialize the scheduler first
  Scheduler* sch = new BasicScheduler;
  sch->init(argc, argv);

  RegisterTask(potrf, void, int, DoubleArray);
  RegisterTask(trsm, void, int, DoubleArray, DoubleArray);
  RegisterTask(syrk, void, int, bool, double, double, DoubleArray, DoubleArray);
  RegisterTask(gemm, void, int, bool, bool, double, double, DoubleArray, DoubleArray, DoubleArray);

  int nblocks = 3;
  int blockSize = 3;
  Matrix A(nblocks, blockSize);

  int ncopies = 1;
  sch->allocateHeap(ncopies);

  for (int i=0; i < ncopies; ++i, sch->nextIter()){
    Task* root = 0;
    if (sch->rank() == 0){
      root = initDag(A);
    }
    sch->run(root);
  }
  sch->deallocateHeap();

  return 0;
}

