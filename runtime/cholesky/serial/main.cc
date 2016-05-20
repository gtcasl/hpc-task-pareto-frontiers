#include <cstdio>
#include <cstring>
#include <vector>
#include <cstdlib>

void
printLowerSymmetric(const char* name, double* A, int n)
{
  printf("%s\n", name);
  double* Aptr = A;
  for (int i=0; i < n; ++i){
    for (int j=0; j <= i; ++j, ++Aptr){
      printf("%10.6f ", *Aptr);
    }
    printf("\n");
  }
  printf("\n");
}

void
printMatrix(const char* name, double* A, int k, int n)
{
  printf("%s\n", name);
  double* Aptr = A;
  for (int i=0; i < k; ++i){
    for (int j=0; j < n; ++j, ++Aptr){
      printf("%10.6f ", *Aptr);
    }
    printf("\n");
  }
  printf("\n");
}


class Matrix
{
 public:
  struct Block {
    int rowStart;
    int colStart;
    int size;
    double* ptr;

    Block(){}

    Block(int r, int c, int s, double* p):
      rowStart(r), colStart(c), size(s), ptr(p)
    {
    }
    
    void print(const char* label)
    {
      printMatrix(label, ptr, size, size);
    }

    void symmetricFill(){
      int rowStop = rowStart + size;
      int colStop = colStart + size;
      int idx = 0;
      double* valptr = ptr;
      for (int i=rowStart; i != rowStop; ++i){
        for (int j=colStart; j != colStop; ++j, ++valptr, ++idx){
          if (i >= j){
            double x = double(0.05 + (idx*7)%4 - (idx*5)%2) / ((idx*11)%5+0.1);
            *valptr = x;
            //printf("Setting point %p to value %f\n", valptr, x);
          } else {
            *valptr = 0;
            //printf("Setting point %p to zero\n", valptr);
          }
        }
      }
    }
  };

  Matrix(int matrixSize, int blockSize) : 
    blockSize_(blockSize), 
    matrixSize_(matrixSize),
    nBlocks_(matrixSize*matrixSize)
  {
    blockStorageSize_ = blockSize * blockSize;

    storageSize_ = nBlocks_*blockStorageSize_;
    storage_ = new double[storageSize_];
    ::memset(storage_, storageSize_*sizeof(double), 0);
    int stride = blockSize*blockSize;
    double* ptr = storage_;
  }

  void print(const char* label = 0){
    if (label) printf("%s\n", label);
    int nrows = blockSize_ * matrixSize_;
    int ncols = nrows;
    int lastRowBlock = 0;
    Block b;
    for (int i=0; i < nrows; ++i){
      int lastColBlock = 0;
      int iBlock = i / blockSize_;
      int iOffset = i % blockSize_;
      if (iBlock != lastRowBlock){
        printf("\n");
        lastRowBlock = iBlock;
      }
      int offset = (iBlock*matrixSize_*blockStorageSize_) + iOffset*blockSize_;
      double* ptr = storage_ + offset;
      for (int j=0; j < ncols; ++j, ++ptr){
        int jBlock = j / blockSize_;
        int jOffset = j % blockSize_;
        if (jBlock != lastColBlock){
          int offset = (iBlock*matrixSize_+jBlock)*blockStorageSize_ + (iOffset*blockSize_ + jOffset);
          printf("    ");
          ptr = storage_ + offset;
          lastColBlock = jBlock;
        }
        //printf("Pointer %p is %f\n", ptr, *ptr);
        printf("%10.6f", *ptr);
      }
      printf("\n");
    }
    printf("\n");
  }

  Block
  block(int i, int j){
    int idx = i*matrixSize_ + j;
    double* ptr = storage_ + idx*blockStorageSize_;
    //printf("Returning block (%d,%d) at offset %d\n", i, j, idx);
    return Block(i*blockSize_, j*blockSize_, blockSize_, ptr);
  }

 private:
  int nBlocks_;
  int blockSize_;
  int matrixSize_;
  double* storage_;
  int storageSize_;
  int blockStorageSize_;
};

extern "C" 
int dgemm(char *transa, char *transb, int *m, int *n, int *k, 
  double *alpha, double *a, int *lda, 
  double *b, int *ldb, double *beta, 
  double *c, int *ldc);


void
multiply(bool ltrans, 
  bool rtrans,
  double alpha,
  double beta,
  Matrix::Block& product, 
  Matrix::Block& left, 
  Matrix::Block& right){
  char lt = ltrans ? 'N' : 'T'; //GD fortran
  char rt = rtrans ? 'N' : 'T'; //GD fortran
  dgemm(&lt, &rt, &product.size, &left.size, &right.size,
    &alpha, left.ptr, &left.size, 
    right.ptr, &right.size, 
    &beta, product.ptr, &product.size);
}


extern "C"
int dsyrk(char *uplo, char *trans, int *n, int *k, 
  double *alpha, double *a, int *lda, double *beta, 
  double *c, int *ldc);

void
symmUpdate(
  bool trans,
  double alpha, 
  double beta,
  Matrix::Block& C,
  Matrix::Block& A)
{
  char ta = trans ? 'N' : 'T';
  char uplo = 'U';
  dsyrk(&uplo, &ta, &A.size, &A.size, 
    &alpha, A.ptr, &A.size,
    &beta, C.ptr, &C.size);

  for (int i=0; i < C.size; ++i){
    for (int j=i+1; j < C.size; ++j){
      int toIdx = i*C.size + j;
      int fromIdx = j*C.size + i;
      C.ptr[toIdx] = C.ptr[fromIdx];
    }
  }
}

extern "C" 
int dpotrf(char *uplo, int *n, double *a, int* lda, int* info);

void
cholesky(Matrix::Block& b)
{
  char uplo = 'U';
  int info;
  dpotrf(&uplo, &b.size, b.ptr, &b.size, &info);
  if (info != 0){
    fprintf(stderr, "FAILURE on DPOTRF: %d\n", info);
    abort();
  }
  double* ptr = b.ptr;
  for (int i=0; i < b.size; ++i){
    for (int j=0; j < b.size; ++j, ++ptr){
      if (j > i) *ptr = 0;
    }
  }
}

extern "C" 
int dtrsm(char *side, char *uplo, char *transa, char *diag, 
  int *m, int *n, 
  double* alpha, double* a, int* lda, 
  double* b, int* ldb);

void
solveX(Matrix::Block& A, Matrix::Block& B)
{
  //solve AX = B
  //B overwrriten with X
  char side = 'R';
  char uplo = 'U';
  char trans = 'N';
  char diag = 'N';
  double alpha = 1.0;
  dtrsm(&side, &uplo, &trans, &diag,
    &B.size, &B.size, &alpha,
    A.ptr, &A.size,
    B.ptr, &B.size);
}

int main(int argc, char** argv)
{
  int nBlocks = 3;
  int nElems = 3;
  char uplo = 'L';

  Matrix A(nBlocks, nElems);
  for (int i=0; i < nBlocks; ++i){
    for (int j=0; j < nBlocks; ++j){
      A.block(i,j).symmetricFill();
    }
  }
  A.print("starting L");

  Matrix B(nBlocks, nElems);
  for (int i=0; i < nBlocks; ++i){
    for (int j=0; j < nBlocks; ++j){
      Matrix::Block Bij = B.block(i,j);
      for (int k=0; k < nBlocks; ++k){
        Matrix::Block Aik = A.block(i,k);
        Matrix::Block Ajk = A.block(j,k);
        multiply(false, true, 1.0, 1.0, Bij, Aik, Ajk);
      }
    }
  }

  B.print("initial");

  for (int k=0; k < nBlocks; ++k){
    Matrix::Block Bkk = B.block(k,k);
    cholesky(Bkk);
    B.print("after diag Cholesky");
    for (int m=k+1; m < nBlocks; ++m){
      Matrix::Block Bmk = B.block(m,k);
      solveX(Bkk, Bmk);
      B.print("eliminate off-diag below");
    }
    for (int n=k+1; n < nBlocks; ++n){
      Matrix::Block Bnn = B.block(n,n);
      Matrix::Block Bnk = B.block(n,k);
      //Bmm -= Bmk * Bmk^T
      //Bnk ends up transposed from the DTRSM solve
      symmUpdate(true, -1.0, 1.0, Bnn, Bnk);
      B.print("after symm update");
      for (int m=n+1; m < nBlocks; ++m){
        Matrix::Block Bmk = B.block(m,k);
        Matrix::Block Bmn = B.block(m,n);
        //because of how awesome BLAS is, Bnk is transposed backwards
        multiply(true, false, -1.0, 1.0, Bmn, Bmk, Bnk);
      }
    }
    B.print("after dgemms");
  }

  B.print("final");
  A.print("starting L");

#if 0
  double* L = new double[n*n];
  double* Lptr = L;
  int idx = 0;
  for (int i=0; i < n; ++i){
    for (int j=0; j < n; ++j, ++idx){
      if (i >= j){
        double x = double(0.05 + (idx*7)%4 - (idx*5)%2) / ((idx*11)%5+0.1);
        printf("A[%d]=x=%f\n", idx, x);
        L[idx] = x;
      } else {
        printf("A[%d]=x=%f\n", idx, 0.0);
        L[idx] = 0;
      }
    }
  }
  printMatrix("initial L", L, n, n);

  double* A = new double[n*n];
  ::memset(A, n*n*sizeof(double), 0);

  double* B = new double[n*n];
  ::memset(B, n*n*sizeof(double), 0);

  char Atrans = 'T';
  char Btrans = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  dgemm(&Atrans, &Btrans, &n, &n, &n, 
    &alpha, L, &n, L, &n, &beta,
    B, &n);
  printMatrix("dgemm test", B, n, n);

  char trans = 'T';
  printMatrix("init A", A, n, n);
  dsyrk(&uplo, &trans, &n, &n, &alpha, L, &n, &beta, A, &n);

  //printLowerSymmetric("dsyrk A", A, n);
  printMatrix("dsyrk A", A, n, n);

  //construct a pos def matrix
  int LDA = n;
  int info;
  dpotrf(&uplo, &n, A, &LDA, &info);
  if (info == 0){
    printf("success\n");
  } else {
    printf("failure %d\n", info);
  }

  printMatrix("final L", A, n, n);
#endif
  return 0;
}

