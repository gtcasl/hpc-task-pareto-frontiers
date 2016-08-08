#ifdef no_mkl
#include <fake_mkl.h>
#else
#include <mkl.h>
#endif
#include <test.h>
#include <matrix.h>

#if CHOLESKY_DEBUG
#define task_debug(...) printf(__VA_ARGS__)
#else
#define task_debug(...) 
#endif

static int static_info;

static enum fxn_id {
  getrf_id,
  gemm_id,
  trsm_left_id, 
  trsm_right_id, 
} test_ids;


void
trsm_left(int k, int m, int size, DoubleArray L, DoubleArray A)
{
  //A overwritten with new U
  task_debug("Solving TRSM, L triangular  A(%d,%d) = L(%d,%d)*U(%d,%d)\n", 
             m, k, m, k, k, k);
  //solve AX = B
  //B overwrriten with X
  char side = 'L';
  char uplo = 'U';
  char trans = 'N';
  char diag = 'N';
  double alpha = 1.0;
  cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
    size, size, alpha,
    L, size,
    A, size);
}

void
trsm_right(int k, int m, int size, DoubleArray A, DoubleArray U)
{
  //A overwritten with new L 
  task_debug("Solving TRSM, U triangular  A(%d,%d) = L(%d,%d)*U(%d,%d)\n", 
             m, k, m, k, k, k);
  //solve AX = B
  //B overwrriten with X
  char side = 'R';
  char uplo = 'U';
  char trans = 'N';
  char diag = 'N';
  double alpha = 1.0;
  cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
    size, size, alpha,
    U, size,
    A, size);
}

void
getrf(int k, int size, IntArray ipiv, DoubleArray A)
{
  task_debug("Running GETRF A(%d,%d)\n", k, k);
  //std::cout << A << std::endl;
  //for(int i = 0; i < size; i++){
  //  std::cout << A[i] << std::endl;
  //}
  char uplo = 'U';
  int info;
  dgetrf(&size, &size, A, &size, ipiv, &info);

  if (info != 0){
    fprintf(stderr, "FAILURE on DGETRF: %d\n", info);
    abort();
  }

}

static void
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
initDag(Matrix& A, IntChunkArray& ipivChunks)
{
  Task* root = 0;

  int nBlocks = A.blockGridSize();
  int blockSize = A.blockSize();

  TaskMap gets(nBlocks);
  TaskMap trs(nBlocks);
  TaskMap gemms(nBlocks);
   
  for (int k=0; k < nBlocks; ++k){
    DoubleArray Akk = A.block(k,k);
    IntArray ipiv = ipivChunks[k];
    Task* diag = new_task(getrf, k, blockSize, ipiv, Akk);
    dep_debug("GETRF (%d,%d) = %p\n", k,k,diag);
    gets(k,k) = diag;
    if (k == 0){
      root = diag;
    } else {
      //the potrf will only directly depend on the most recent syrk
      diag->dependsOn(gemms(k,k));
    }
    for (int m=k+1; m < nBlocks; ++m){
      DoubleArray Amk = A.block(m,k);
      Task* solve = new_task(trsm_left, k, m, blockSize, Akk, Amk);
      dep_debug(" TRSM (%d,%d) = %p\n", m, k, solve);
      trs(m,k) = solve;
      solve->dependsOn(gets(k,k));
      if (m >= 2){
        solve->dependsOn(gemms(m,k));
      }

      DoubleArray Akm = A.block(k,m);
      solve = new_task(trsm_right, k, m, blockSize, Akk, Akm);
      dep_debug(" TRSM (%d,%d) = %p\n", k, m, solve);
      trs(k,m) = solve;
      solve->dependsOn(gets(k,k));
      if (m >= 2){
        solve->dependsOn(gemms(k,m));
      }
    }
    for (int m=k+1; m < nBlocks; ++m){
      DoubleArray Amk = A.block(m,k);
      for (int n=k+1; n < nBlocks; ++n){
        DoubleArray Akn = A.block(k,n);
        DoubleArray Amn = A.block(m,n);
        Task* mult = new_task(gemm, m, n, k, blockSize, true, false, -1.0, 1.0, Amn, Amk, Akn);
        dep_debug(" GEMM (%d,%d) = %p\n", m,n,mult);
        gemms(m,n) = mult;
        mult->dependsOn(trs(m,k));
        mult->dependsOn(trs(k,n));
      }
    }
  }
  return root;
}

int lu(int argc, char** argv)
{
  //ALWAYS Initialize the scheduler first
#ifdef basic_scheduler
  Scheduler* sch = new BasicScheduler;
#else
  Scheduler* sch = new AdvancedScheduler;
#endif
  sch->init(argc, argv);
  // disable dynamic thread adjustment in MKL
#ifndef no_mkl
  mkl_set_dynamic(0);
#endif

  RegisterTask(getrf, void, int,
    int, IntArray, DoubleArray);
  RegisterTask(trsm_left,  void, int, int,
    int, DoubleArray, DoubleArray);
  RegisterTask(trsm_right,  void, int, int,
    int, DoubleArray, DoubleArray);
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
  IntArray ipiv(blockSize*nBlocks);
  IntChunkArray ipivChunks(nBlocks);



  int ncopies = 1;
  sch->allocateHeap(ncopies);

  for (int i=0; i < nBlocks; ++i){
    ipivChunks[i] = ipiv.offset(i*blockSize); 
  }

  for (int c=0; c < ncopies; ++c, sch->nextIter()){
    if (sch->rank() == 0){

      for (int i=0; i < nBlocks; ++i){
        for (int j=0; j < nBlocks; ++j){
          A.symmetricFill(i,j);
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  sch->resetIter();

  for (int iter=0; iter < ncopies; ++iter, sch->nextIter()){
    Task* root = 0;
    if (sch->rank() == 0){
      root = initDag(A, ipivChunks);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    sch->run(root);
  }

  fflush(stdout);

  //A.print("factorized A");


  sch->deallocateHeap();
  sch->finalize();

  return 0;
}

RegisterTest("lu", lu);
