#ifdef no_mkl
#include <fake_mkl.h>
#else
#include <mkl.h>
#endif
#include <test.h>

typedef Buffer<double> DoublePtr;
typedef Buffer<double> DoubleArray;
typedef Buffer<int> IntArray;
typedef Buffer<IntArray> IntChunkArray;
typedef Buffer<DoubleArray> DoubleChunkArray;


#define debug(x) printf(#x "\n");
#undef debug
#define debug(x) 

struct config
{
  int nchunks;
  int nnz_per_row;
  int nrows_per_chunk;
  int nrows;
  int ncols;
  int niter;
  void print(){
    std::cout << "nchunks: " << nchunks << "\n"
              << "nnz_per_row: " << nnz_per_row << "\n"
              << "nrows: " << nrows<< "\n"
              << "ncols: " << ncols<< "\n"
              << "niter: " << niter<< "\n";
  }
};

static enum fxn_id {
  myddot_id,
  start_id,
  mydaxpy_id,
  dxapy_id,
  spmv_id,
  copy_id,
  subtract_id,
  sum_contribs_id,
  comp_alpha_id,
  comp_beta_id,
  assign_id
} test_ids;

extern void
generate_problem_27pt(
  int nx, int ny, int nz,
  int chunkSize,
  DoubleArray A,
  IntArray nonzeros,
  IntArray nnzPerRow,
  DoubleChunkArray bChunks,
  DoubleChunkArray xChunks,
  DoubleChunkArray AChunks,
  IntChunkArray    nonzerosChunks
);

extern void
multiply(
  int nrows,
  DoubleArray vector,
  DoubleArray residual,
  DoubleArray matrix,
  IntArray nnzPerRow,
  IntArray nonzerosInRow);

void sum_contribs(int n, DoubleArray contribs, DoublePtr result)
{
  debug(start sum_contribs);
  for(int i = 0; i < n; ++i){
    *result += contribs[i];
  }
  debug(end sum_contribs);
}

//Compute Ax = y
void spmv(int nrows, DoubleArray A, DoubleArray x, DoubleArray y, IntArray nnzPerRow, IntArray nonzerosInRow)
{
  debug(start spmv);
  //cblas_dspmv(CblasRowMajor, CblasUpper, nrows, 1.0, A, x, 1, 0.0, y, 1);
  multiply(nrows, x, y, A, nnzPerRow, nonzerosInRow);
  debug(end spmv);
}


/**
 * Compute C = A - B
 */
static void subtract(int n, DoubleArray C, DoubleArray A, DoubleArray B)
{
  debug(start subtract);
  vdSub(n, A, B, C);
  debug(end subtract);
}

static void myddot(int n, DoubleArray A, DoubleArray B, DoublePtr result)
{
  debug(start myddot);
  *result = cblas_ddot(n, A, 1, B, 1);
  debug(end myddot);
}

/**
 * @brief daxpy  Computes Y = A*X + Y
 * @param n
 * @param a
 * @param x
 * @param y
 */
void mydaxpy(int n, double a, DoubleArray x, DoubleArray y)
{
  debug(start mydaxpy);
  cblas_daxpy(n, a, x, 1, y, 1);
  debug(end mydaxpy);
}

/**
 * @brief dxapy Computes Y = X + A*Y
 * @param n
 * @param a
 * @param x
 * @param y
 */
void dxapy(int n, double a, DoubleArray x, DoubleArray y)
{
  debug(start dxapy);
  //double scale = (1+a);
  //scale vector Y by (1+a) - then add in 1.0 times x
  //cblas_dscal(n, a, y, 1);
  //cblas_daxpy(n, 1.0, x, 1, y, 1);
  cblas_daxpby(n, 1.0, x, 1, a, y, 1);
  debug(end dxapy);
}

void copy(int n, DoubleArray dst, DoubleArray src)
{
  debug(start copy);
  cblas_dcopy(n, src, 1, dst, 1);
  debug(end copy);
}

void comp_alpha(double rsq, double pAp, DoublePtr result)
{
  debug(start comp_alpha);
  *result = rsq / pAp;
  debug(end comp_alpha);
}

void comp_beta(double rsq, double rsqNextIter, DoublePtr result)
{
  debug(start comp_beta);
  *result = rsqNextIter / rsq;
  debug(end comp_beta);
}

void assign(DoublePtr dst, DoublePtr src)
{
  *dst = *src;
}


void start(int ignore)
{
  debug(start start);
  //no op for organizational purposes
  debug(end start);
}

Task*
initDag(config cfg,
  DoubleArray A,
  DoubleArray Ap,
  DoubleArray p,
  DoubleArray r,
  DoubleArray x,
  DoubleArray pApContribs,
  DoubleArray RsqContribs,
  DoublePtr alpha,
  DoublePtr beta,
  DoublePtr pAp,
  DoublePtr Rsq,
  DoublePtr RsqNextIter,
  DoubleChunkArray  AChunks,
  DoubleChunkArray  ApChunks,
  DoubleChunkArray  bChunks,
  DoubleChunkArray  pChunks,
  DoubleChunkArray  rChunks,
  DoubleChunkArray  xChunks,
  IntChunkArray  nnzPerRowChunks,
  IntChunkArray  nonzerosChunks
)
{
  Task* root = new_task(start, 1);
  Task* sumRsq  = new_task(sum_contribs, cfg.nchunks, RsqContribs, Rsq);
  std::vector<Task*> lastWavefront(cfg.nchunks);
  for (int i=0; i < cfg.nchunks; ++i){
    //have to start by computing Ax
    Task* compAx  = new_task(spmv, cfg.nrows_per_chunk, AChunks[i], x, ApChunks[i], nnzPerRowChunks[i], nonzerosChunks[i]);
    Task* compR0  = new_task(subtract, cfg.nrows_per_chunk, rChunks[i], bChunks[i], ApChunks[i]);
    Task* compP0  = new_task(copy, cfg.nrows_per_chunk, pChunks[i], rChunks[i]);
    Task* compRsq = new_task(myddot,  cfg.nrows_per_chunk, rChunks[i], rChunks[i], RsqContribs.offset(i));

    compAx->dependsOn(root);
    compP0->dependsOn(compR0);
    compR0->dependsOn(compAx);
    compRsq->dependsOn(compR0);
    sumRsq->dependsOn(compRsq);
    lastWavefront[i] = compP0;
  }

  Task* iterPrecursor = sumRsq;
  //now we enter the loop
  for (int iter=0; iter < cfg.niter; ++iter){
    Task* sum_pAp = new_task(sum_contribs, cfg.nchunks, pApContribs, pAp);
    for (int i=0; i < cfg.nchunks; ++i){
      Task* compAp    = new_task(spmv,  cfg.nrows_per_chunk, AChunks[i], p, ApChunks[i], nnzPerRowChunks[i], nonzerosChunks[i]);
      Task* comp_pAp  = new_task(myddot,  cfg.nrows_per_chunk, ApChunks[i], pChunks[i], pApContribs.offset(i));
      compAp->dependsOn(lastWavefront[i]);
      comp_pAp->dependsOn(compAp);
      sum_pAp->dependsOn(comp_pAp);
    }

    Task* compAlpha = new_task(comp_alpha, *alpha, *pAp, Rsq);
    compAlpha->dependsOn(iterPrecursor);
    compAlpha->dependsOn(sum_pAp);

    sumRsq = new_task(sum_contribs, cfg.nchunks, RsqContribs, RsqNextIter);
    for (int i=0; i < cfg.nchunks; ++i){
      Task* compX     = new_task(mydaxpy, cfg.nrows_per_chunk, *alpha, pChunks[i], xChunks[i]);
      Task* compR     = new_task(mydaxpy, cfg.nrows_per_chunk, -(*alpha), ApChunks[i], rChunks[i]);
      Task* compRsq   = new_task(myddot,  cfg.nrows_per_chunk, rChunks[i], rChunks[i], RsqContribs.offset(i));
      compX->dependsOn(compAlpha);
      compR->dependsOn(compAlpha);
      compRsq->dependsOn(compR);
      sumRsq->dependsOn(compRsq);
    }

    Task* compBeta  = new_task(comp_beta, *Rsq, *RsqNextIter, beta);
    Task* setEqual  = new_task(assign, RsqNextIter, Rsq);
    compBeta->dependsOn(sumRsq);
    setEqual->dependsOn(sumRsq);

    for (int i=0; i < cfg.nchunks; ++i){
      Task* compP  = new_task(dxapy, cfg.nrows_per_chunk, *beta, rChunks[i], pChunks[i]);
      compP->dependsOn(compBeta);
      lastWavefront[i] = compP;
    }

    iterPrecursor = setEqual;
  }

  return root;
}


int cg(int argc, char** argv)
{
  //ALWAYS Initialize the scheduler first
  Scheduler* sch = new BasicScheduler;
  sch->init(argc, argv);

  // disable dynamic thread adjustment in MKL
  mkl_set_dynamic(0);

  RegisterTask(start,void,int);
  RegisterTask(sum_contribs, void, int, DoubleArray, DoublePtr);
  RegisterTask(spmv, void, int, DoubleArray, DoubleArray, DoubleArray, IntArray, IntArray);
  RegisterTask(myddot, void, int, DoubleArray, DoubleArray, DoublePtr);
  RegisterTask(mydaxpy, void, int, double, DoubleArray, DoubleArray);
  RegisterTask(dxapy, void, int, double, DoubleArray, DoubleArray);
  RegisterTask(comp_alpha, void, double, double, DoublePtr);
  RegisterTask(comp_beta, void, double, double, DoublePtr);
  RegisterTask(assign, void, DoublePtr, DoublePtr);
  RegisterTask(copy,void,int,DoubleArray,DoubleArray);
  RegisterTask(subtract,void,int,DoubleArray,DoubleArray,DoubleArray);

  config cfg;
  if(argc != 6){
    if (sch->rank() == 0){
      std::cerr << "Usage: " << argv[1] << " <nx> <ny> <nz> <nchunks>" << std::endl;
    }
    return -1;
  }

  int nx = atoi(argv[2]), ny = atoi(argv[3]), nz = atoi(argv[4]);
  int ncopies = 1;
  int nrows = nx*ny*nz;
  int nnz_per_row = 27;
  int nchunks = atoi(argv[5]);
  cfg.nchunks = nchunks;
  cfg.nnz_per_row = nrows;
  cfg.nrows = nx*ny*nz;
  cfg.ncols = cfg.nrows;
  cfg.niter = 1;
  cfg.nrows_per_chunk = nrows / nchunks;

  if(sch->rank() == 0){
    cfg.print();
  }

  int chunkSize = nrows / nchunks;
  if (nrows % nchunks){
    std::cerr << "number of chunks does not divide nx*ny*nz - reconfigure" << std::endl;
    return 1;
  }

  DoubleArray p(nrows);
  DoubleArray x(nrows);
  DoubleArray r(nrows);
  DoubleArray A(nrows*nnz_per_row);
  DoubleArray Ap(nrows);
  IntArray    nnzPerRow(nrows);
  IntArray    nonzeros(nrows*nnz_per_row);
  DoubleArray RsqContribs(nchunks);
  DoubleArray pApContribs(nchunks);
  DoublePtr   Alpha(1);
  DoublePtr   Beta(1);
  DoublePtr   pAp(1);
  DoublePtr   Rsq(1);
  DoublePtr   RsqNextIter(1);

  DoubleChunkArray pChunks(nchunks);
  DoubleChunkArray bChunks(nchunks);
  DoubleChunkArray xChunks(nchunks);
  DoubleChunkArray rChunks(nchunks);
  DoubleChunkArray AChunks(nchunks);
  DoubleChunkArray ApChunks(nchunks);
  IntChunkArray    nnzPerRowChunks(nchunks);
  IntChunkArray    nonzerosChunks(nchunks);

  sch->allocateHeap(ncopies);

  for (int iters=0; iters < ncopies; ++iters, sch->nextIter()){
    for (int i=0; i < nchunks; ++i){
      pChunks[i] = p.offset(i*chunkSize);
      xChunks[i] = x.offset(i*chunkSize);
      rChunks[i] = r.offset(i*chunkSize);
      nnzPerRowChunks[i] = nnzPerRow.offset(i*chunkSize);
    }
    generate_problem_27pt(nx,ny,nz,chunkSize,A,nonzeros,nnzPerRow,bChunks,xChunks,AChunks,nonzerosChunks);
  }

  for (int i=0; i < ncopies; ++i, sch->nextIter()){
    Task* root = 0;
    if (sch->rank() == 0){
      root = initDag(cfg,
        A, Ap, p, r, x, //full arrays
        pApContribs, RsqContribs, //reduction arrays
        Alpha, Beta, pAp, Rsq, RsqNextIter, //double pointers
        AChunks, ApChunks, bChunks, pChunks, rChunks, xChunks, //double chunk arrays,
        nnzPerRowChunks, nonzerosChunks //int chunk arrays
      );
    }
    sch->run(root);
    sch->stop();
  }
  sch->deallocateHeap();
  sch->finalize();

  return 0;
}

RegisterTest("cg", cg);

