#ifndef _cholesky_h_
#define _cholesky_h_

enum fxn_id {
  potrf_id,
  gemm_id,
  syrk_id,
  trsm_id,
  dopotrf_id,
  dogemm_id,
  dosyrk_id,
  dotrsm_id,
  potrfprofile_id,
  syrkprofile_id,
  trsmprofile_id,
  gemmprofile_id,
} ;

void trsm(int k, int m, int size, int lda, double* A, double* B);
void potrf(int k, int size, int lda, double* A);
void syrk(int k, int n, int size, int lda, double* C, double* A);
void gemm(int m, int n, int k, int size, int lda, bool ltrans,  bool rtrans, double alpha, double beta, double* product,  double* left,  double* right);

class Matrix
{
 private:
  int blockOffset(int row, int col){
    int lda = blockGridSize_ * blockSize_;
    return lda * blockSize_ * col + blockSize_ * row;
  }

 public:
  void
  randomFill(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0,1);
    int nrows = blockSize_ * blockGridSize_;
//#pragma omp parallel for
    for(int i = 0; i < nrows * nrows; i++){
      //double x = double(0.05 + (i*7)%4 - (i*5)%2) / ((i*11)%5+0.1);
      double x = dis(gen);
      storage_[i] = x;
    }
  }

  void 
  symmetricFill(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0,1);
    int nrows = blockSize_ * blockGridSize_;
//#pragma omp parallel for
    for(int i = 0; i < nrows; i++){
      for(int j = 0; j < nrows; j++){
        int idx = i*nrows + j;
        if(i == j){
          double x = dis(gen);
          storage_[idx] = x;
        } else {
          storage_[idx] = 0.0;
        }
      }
    }
  }

  Matrix(int matrixSize, int blockSize) : 
    blockSize_(blockSize), 
    blockGridSize_(matrixSize)
  {
    storage_ = (double*) mkl_malloc(matrixSize * matrixSize * blockSize * blockSize * sizeof(double), 64);
  }

  void print(const char* label = 0){
    if (label) printf("%s\n", label);
    int nrows = blockSize_ * blockGridSize_;
    for(int i = 0; i < nrows; i++){
      if(i > 0 && i % blockSize_ == 0){
        printf("\n");
      }
      for(int j = 0; j < nrows; j++){
        if(j > 0 && j % blockSize_ == 0){
          printf("   ");
        }
        auto val = storage_[j*nrows + i];
        printf("%6.4f ", val);
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

  double*
  block(int i, int j){
    return storage_ + blockOffset(i,j);
  }

  double*
  storageAddr(){
    return storage_;
  }

  int blockSize_;
  int blockGridSize_;
  double* storage_;
};


#endif
