#ifndef _matrix_h_
#define _matrix_h_

class Matrix
{
 private:
  int blockOffset(int row, int col){
    int lda = blockGridSize_ * blockSize_;
    return lda * blockSize_ * col + blockSize_ * row;
  }

 public:
  void randomFill();

  void symmetricFill();

  Matrix(int matrixSize, int blockSize);

  void print(const char* label = 0);

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
