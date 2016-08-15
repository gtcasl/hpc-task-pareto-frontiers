#ifdef no_mkl
#include <fake_mkl.h>
#else
#include <mkl.h>
#endif
#include <test.h>

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

enum fxn_id {
  dot_id,
  start_id,
  axpy_id,
  axpyprofile_id,
  xapy_id,
  xapyprofile_id,
  spmv_id,
  copy_id,
  subtract_id,
  sum_contribs_id,
  comp_alpha_id,
  comp_beta_id,
  assign_id
};

extern void
generate_problem_27pt(
  int nx, int ny, int nz,
  int chunkSize,
  double* A,
  int* nonzeros,
  int* nnzPerRow,
  double** AChunks,
  int**    nonzerosChunks
);

void sum_contribs(int n, double* contribs, double* result)
{
  debug(start sum_contribs);
#pragma offload target(mic:0) \
  in(contribs:length(0) alloc_if(0) free_if(0)) \
  in(result:length(0) alloc_if(0) free_if(0))
{
  double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
  for(int i = 0; i < n; ++i){
    sum += contribs[i];
  }
  *result = sum;
}
  debug(end sum_contribs);
}

//Compute Ax = y
void spmv(int nrows, double* A, double* x, double* y, int* nnzPerRow, int* nonzerosInRow)
{
  debug(start spmv);
#pragma offload target(mic:0) \
  in(y : length(0) alloc_if(0) free_if(0)) \
  in(x : length(0) alloc_if(0) free_if(0)) \
  in(nnzPerRow : length(0) alloc_if(0) free_if(0)) \
  in(nonzerosInRow : length(0) alloc_if(0) free_if(0)) \
  in(A : length(0) alloc_if(0) free_if(0)) 
{
  double* v = x;
  double* r = y;
  int* nonzeros = nonzerosInRow;
  #pragma omp parallel for
  for (int row=0; row < nrows; ++row){
    double res = 0; 
    int nnz = nnzPerRow[row];
    for (int i=0; i < nnz; ++i){
      int col = nonzeros[i];
      res += A[i] * v[col];
    }
    A += nnz;
    nonzeros += nnz;
    r[row] = res;
  }
}
  debug(end spmv);
}


/**
 * Compute C = A - B
 */
void subtract(int n, double* C, double* A, double* B)
{
  if(!C){std::cout << "C" << std::endl;}
  if(!A){std::cout << "A" << std::endl;}
  if(!B){std::cout << "B" << std::endl;}
  debug(start subtract);
#pragma offload target(mic:0) \
  in(A:length(0) alloc_if(0) free_if(0)) \
  in(B:length(0) alloc_if(0) free_if(0)) \
  in(C:length(0) alloc_if(0) free_if(0)) 
{
  vdSub(n, A, B, C);
}

  debug(end subtract);
}

void dot(int n, double* A, double* B, double* result)
{
  debug(start dot);
#pragma offload target(mic:0) \
  in(A:length(0) alloc_if(0) free_if(0)) \
  in(B:length(0) alloc_if(0) free_if(0)) \
  in(result:length(0) alloc_if(0) free_if(0))
{
  double res;
  res = cblas_ddot(n, A, 1, B, 1);
  *result = res;
}
  debug(end dot);
}

/**
 * @brief daxpy  Computes Y = A*X + Y
 * @param n
 * @param a
 * @param x
 * @param y
 */
void axpyprofile(int n, double a, double* x, double* y)
{
  debug(start axpy);
#pragma offload target(mic:0) \
  in(x:length(0) alloc_if(0) free_if(0)) \
  in(y:length(0) alloc_if(0) free_if(0))
{
  double* ycopy = (double*) mkl_malloc(n * sizeof(double), 64);
#pragma omp parallel for
  for(int i = 0; i < n; i++){
    ycopy[i] = y[i];
  }
  cblas_daxpy(n, a, x, 1, ycopy, 1);
  mkl_free(ycopy);
}
  debug(end axpy);
}

/**
 * @brief daxpy  Computes Y = A*X + Y
 * @param n
 * @param a
 * @param x
 * @param y
 */
void axpy(int n, double a, double* x, double* y)
{
  debug(start axpy);
#pragma offload target(mic:0) \
  in(x:length(0) alloc_if(0) free_if(0)) \
  in(y:length(0) alloc_if(0) free_if(0))
{
  cblas_daxpy(n, a, x, 1, y, 1);
}
  debug(end axpy);
}

/**
 * @brief dxapy Computes Y = X + A*Y
 * @param n
 * @param a
 * @param x
 * @param y
 */
void xapyprofile(int n, double a, double* x, double* y)
{
  debug(start dxapy);
  //double scale = (1+a);
  //scale vector Y by (1+a) - then add in 1.0 times x
  //cblas_dscal(n, a, y, 1);
  //cblas_daxpy(n, 1.0, x, 1, y, 1);
#pragma offload target(mic:0) \
  in(x:length(0) alloc_if(0) free_if(0)) \
  in(y:length(0) alloc_if(0) free_if(0))
{
  double* ycopy = (double*) mkl_malloc(n * sizeof(double), 64);
#pragma omp parallel for
  for(int i = 0; i < n; i++){
    ycopy[i] = y[i];
  }
  cblas_daxpby(n, 1.0, x, 1, a, y, 1);
  mkl_free(ycopy);
}
  debug(end dxapy);
}

/**
 * @brief dxapy Computes Y = X + A*Y
 * @param n
 * @param a
 * @param x
 * @param y
 */
void xapy(int n, double a, double* x, double* y)
{
  debug(start dxapy);
  //double scale = (1+a);
  //scale vector Y by (1+a) - then add in 1.0 times x
  //cblas_dscal(n, a, y, 1);
  //cblas_daxpy(n, 1.0, x, 1, y, 1);
#pragma offload target(mic:0) \
  in(x:length(0) alloc_if(0) free_if(0)) \
  in(y:length(0) alloc_if(0) free_if(0))
{
  cblas_daxpby(n, 1.0, x, 1, a, y, 1);
}
  debug(end dxapy);
}

void copy(int n, double* dst, double* src)
{
  debug(start copy);
#pragma offload target(mic:0) \
  in(dst:length(0) alloc_if(0) free_if(0)) \
  in(src:length(0) alloc_if(0) free_if(0))
{
  cblas_dcopy(n, src, 1, dst, 1);
}
  debug(end copy);
}

void comp_alpha(double* rsq, double* pAp, double* result)
{
  debug(start comp_alpha);
#pragma offload target(mic:0) \
  in(rsq:length(0) alloc_if(0) free_if(0)) \
  in(pAp:length(0) alloc_if(0) free_if(0)) \
  in(result:length(0) alloc_if(0) free_if(0))
{
  *result = *rsq / *pAp;
}
  debug(end comp_alpha);
}

void comp_beta(double* rsq, double* rsqNextIter, double* result)
{
  debug(start comp_beta);
#pragma offload target(mic:0) \
  in(rsq:length(0) alloc_if(0) free_if(0)) \
  in(rsqNextIter:length(0) alloc_if(0) free_if(0)) \
  in(result:length(0) alloc_if(0) free_if(0))
{
  *result = *rsqNextIter / *rsq;
}
  debug(end comp_beta);
}

void assign(double* dst, double* src)
{
#pragma offload target(mic:0) \
  in(src:length(0) alloc_if(0) free_if(0)) \
  in(dst:length(0) alloc_if(0) free_if(0))
{
  *dst = *src;
}
}


void start(int)
{
  debug(start start);
  //no op for organizational purposes
  debug(end start);
}

Task*
initProfilingDag(config cfg,
  double* p,
  double* x,
  double* pApContribs,
  double* RsqContribs,
  double* alpha,
  double* beta,
  double* pAp,
  double* Rsq,
  double* RsqNextIter,
  double**  AChunks,
  double**  ApChunks,
  double**  bChunks,
  double**  pChunks,
  double**  rChunks,
  double**  xChunks,
  int**  nnzPerRowChunks,
  int**  nonzerosChunks
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
    Task* compRsq = new_task(dot,  cfg.nrows_per_chunk, rChunks[i], rChunks[i], RsqContribs + i);

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
      Task* comp_pAp  = new_task(dot,  cfg.nrows_per_chunk, ApChunks[i], pChunks[i], pApContribs + i);
      compAp->dependsOn(lastWavefront[i]);
      comp_pAp->dependsOn(compAp);
      sum_pAp->dependsOn(comp_pAp);
    }

    Task* compAlpha = new_task(comp_alpha, pAp, Rsq, alpha);
    compAlpha->dependsOn(iterPrecursor);
    compAlpha->dependsOn(sum_pAp);

    sumRsq = new_task(sum_contribs, cfg.nchunks, RsqContribs, RsqNextIter);
    for (int i=0; i < cfg.nchunks; ++i){
      Task* compXprofile     = new_task(axpyprofile, cfg.nrows_per_chunk, *alpha, pChunks[i], xChunks[i]);
      Task* compX     = new_mutating_task(axpy, cfg.nrows_per_chunk, *alpha, pChunks[i], xChunks[i]);
      compX->dependsOn(compXprofile);
      Task* compRprofile     = new_task(axpyprofile, cfg.nrows_per_chunk, -(*alpha), ApChunks[i], rChunks[i]);
      Task* compR     = new_mutating_task(axpy, cfg.nrows_per_chunk, -(*alpha), ApChunks[i], rChunks[i]);
      compR->dependsOn(compRprofile);
      Task* compRsq   = new_task(dot,  cfg.nrows_per_chunk, rChunks[i], rChunks[i], RsqContribs + i);
      compXprofile->dependsOn(compAlpha);
      compRprofile->dependsOn(compAlpha);
      compRsq->dependsOn(compR);
      sumRsq->dependsOn(compRsq);
    }

    Task* compBeta  = new_task(comp_beta, Rsq, RsqNextIter, beta);
    Task* setEqual  = new_task(assign, RsqNextIter, Rsq);
    compBeta->dependsOn(sumRsq);
    setEqual->dependsOn(sumRsq);

    for (int i=0; i < cfg.nchunks; ++i){
      Task* compPprofile  = new_task(xapyprofile, cfg.nrows_per_chunk, *beta, rChunks[i], pChunks[i]);
      Task* compP  = new_mutating_task(xapy, cfg.nrows_per_chunk, *beta, rChunks[i], pChunks[i]);
      compP->dependsOn(compPprofile);
      compPprofile->dependsOn(compBeta);
      lastWavefront[i] = compP;
    }

    iterPrecursor = setEqual;
  }

  return root;
}

Task*
initDag(config cfg,
  double* p,
  double* x,
  double* pApContribs,
  double* RsqContribs,
  double* alpha,
  double* beta,
  double* pAp,
  double* Rsq,
  double* RsqNextIter,
  double**  AChunks,
  double**  ApChunks,
  double**  bChunks,
  double**  pChunks,
  double**  rChunks,
  double**  xChunks,
  int**  nnzPerRowChunks,
  int**  nonzerosChunks
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
    Task* compRsq = new_task(dot,  cfg.nrows_per_chunk, rChunks[i], rChunks[i], RsqContribs + i);

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
      Task* comp_pAp  = new_task(dot,  cfg.nrows_per_chunk, ApChunks[i], pChunks[i], pApContribs + i);
      compAp->dependsOn(lastWavefront[i]);
      comp_pAp->dependsOn(compAp);
      sum_pAp->dependsOn(comp_pAp);
    }

    Task* compAlpha = new_task(comp_alpha, pAp, Rsq, alpha);
    compAlpha->dependsOn(iterPrecursor);
    compAlpha->dependsOn(sum_pAp);

    sumRsq = new_task(sum_contribs, cfg.nchunks, RsqContribs, RsqNextIter);
    for (int i=0; i < cfg.nchunks; ++i){
      Task* compX     = new_task(axpy, cfg.nrows_per_chunk, *alpha, pChunks[i], xChunks[i]);
      Task* compR     = new_task(axpy, cfg.nrows_per_chunk, -(*alpha), ApChunks[i], rChunks[i]);
      Task* compRsq   = new_task(dot,  cfg.nrows_per_chunk, rChunks[i], rChunks[i], RsqContribs + i);
      compX->dependsOn(compAlpha);
      compR->dependsOn(compAlpha);
      compRsq->dependsOn(compR);
      sumRsq->dependsOn(compRsq);
    }

    Task* compBeta  = new_task(comp_beta, Rsq, RsqNextIter, beta);
    Task* setEqual  = new_task(assign, RsqNextIter, Rsq);
    compBeta->dependsOn(sumRsq);
    setEqual->dependsOn(sumRsq);

    for (int i=0; i < cfg.nchunks; ++i){
      Task* compP  = new_task(xapy, cfg.nrows_per_chunk, *beta, rChunks[i], pChunks[i]);
      compP->dependsOn(compBeta);
      lastWavefront[i] = compP;
    }

    iterPrecursor = setEqual;
  }

  return root;
}

int cgprofiling(Scheduler* sch, int argc, char** argv){
  config cfg;
  if(argc != 5){
    std::cerr << "Usage: " << argv[0] << " <nx> <ny> <nz> <nchunks>" << std::endl;
    return -1;
  }

  int nx = atoi(argv[1]), ny = atoi(argv[2]), nz = atoi(argv[3]);
  int nrows = nx*ny*nz;
  int nnz_per_row = 27;
  int nchunks = atoi(argv[4]);
  cfg.nchunks = nchunks;
  cfg.nnz_per_row = nnz_per_row;
  cfg.nrows = nrows;
  cfg.ncols = cfg.nrows;
  cfg.niter = 1;
  cfg.nrows_per_chunk = nrows / nchunks;

  cfg.print();

  int chunkSize = nrows / nchunks;
  if (nrows % nchunks){
    std::cerr << "number of chunks does not divide nx*ny*nz - reconfigure" << std::endl;
    return 1;
  }

  double* p = (double*) mkl_malloc(nrows * sizeof(double), 64);
  double* x = (double*) mkl_malloc(nrows * sizeof(double), 64);
  double* r = (double*) mkl_malloc(nrows * sizeof(double), 64);
  double* A = (double*) mkl_malloc(nrows*nnz_per_row * sizeof(double), 64);
  double* Ap = (double*) mkl_malloc(nrows * sizeof(double), 64);
  int*    nnzPerRow = (int*) mkl_malloc(nrows * sizeof(int), 64);
  int*    nonzeros = (int*) mkl_malloc(nrows*nnz_per_row * sizeof(int), 64);
  double* RsqContribs = (double*) mkl_malloc(nchunks * sizeof(double), 64);
  double* pApContribs = (double*) mkl_malloc(nchunks * sizeof(double), 64);
  double* Alpha = (double*) mkl_malloc(sizeof(double), 64);
  double* Beta = (double*) mkl_malloc(sizeof(double), 64);
  double* pAp = (double*) mkl_malloc(sizeof(double), 64);
  double* Rsq = (double*) mkl_malloc(sizeof(double), 64);
  double* RsqNextIter = (double*) mkl_malloc(sizeof(double), 64);
  double* b = (double*) mkl_malloc(nrows * sizeof(double), 64);
  for(int i = 0; i < nrows; i++){
    b[i] = i;
  }

  double** pChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);
  double** bChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);
  double** xChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);
  double** rChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);

  double** AChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);
  double** ApChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);
  int**    nnzPerRowChunks = (int**) mkl_malloc(nchunks * sizeof(int*), 64);
  int**    nonzerosChunks = (int**) mkl_malloc(nchunks * sizeof(int*), 64);

  for (int i=0; i < nchunks; ++i){
    pChunks[i] = p + i*chunkSize;
    xChunks[i] = x + i*chunkSize;
    rChunks[i] = r + i*chunkSize;
    bChunks[i] = b + i*chunkSize;
    nnzPerRowChunks[i] = nnzPerRow + i*chunkSize;
    ApChunks[i] = Ap + i*chunkSize;
  }
  generate_problem_27pt(nx,ny,nz,chunkSize,A,nonzeros,nnzPerRow,AChunks,nonzerosChunks);

#pragma offload_transfer target(mic:0) \
  in(p : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(x : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(r : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(b : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(A : length(nrows*nnz_per_row) align(64) alloc_if(1) free_if(0)) \
  in(Ap : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(nnzPerRow : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(nonzeros : length(nrows*nnz_per_row) align(64) alloc_if(1) free_if(0)) \
  in(RsqContribs : length(nchunks) align(64) alloc_if(1) free_if(0)) \
  in(pApContribs : length(nchunks) align(64) alloc_if(1) free_if(0)) \
  in(Alpha : length(1) align(64) alloc_if(1) free_if(0)) \
  in(Beta : length(1) align(64) alloc_if(1) free_if(0)) \
  in(pAp : length(1) align(64) alloc_if(1) free_if(0)) \
  in(Rsq : length(1) align(64) alloc_if(1) free_if(0)) \
  in(RsqNextIter : length(1) align(64) alloc_if(1) free_if(0))

  Task* root = 0;
  root = initProfilingDag(cfg,
                 p, x, //full arrays
                 pApContribs, RsqContribs, //reduction arrays
                 Alpha, Beta, pAp, Rsq, RsqNextIter, //double pointers
                 AChunks, ApChunks, bChunks, pChunks, rChunks, xChunks, //double chunk arrays,
                 nnzPerRowChunks, nonzerosChunks //int chunk arrays
                );
  sch->run(root);

  return 0;
}

int cg(Scheduler* sch, int argc, char** argv){
  config cfg;
  if(argc != 5){
    std::cerr << "Usage: " << argv[0] << " <nx> <ny> <nz> <nchunks>" << std::endl;
    return -1;
  }

  int nx = atoi(argv[1]), ny = atoi(argv[2]), nz = atoi(argv[3]);
  int nrows = nx*ny*nz;
  int nnz_per_row = 27;
  int nchunks = atoi(argv[4]);
  cfg.nchunks = nchunks;
  cfg.nnz_per_row = nnz_per_row;
  cfg.nrows = nrows;
  cfg.ncols = cfg.nrows;
  cfg.niter = 1;
  cfg.nrows_per_chunk = nrows / nchunks;

  cfg.print();

  int chunkSize = nrows / nchunks;
  if (nrows % nchunks){
    std::cerr << "number of chunks does not divide nx*ny*nz - reconfigure" << std::endl;
    return 1;
  }

  double* p = (double*) mkl_malloc(nrows * sizeof(double), 64);
  double* x = (double*) mkl_malloc(nrows * sizeof(double), 64);
  double* r = (double*) mkl_malloc(nrows * sizeof(double), 64);
  double* A = (double*) mkl_malloc(nrows*nnz_per_row * sizeof(double), 64);
  double* Ap = (double*) mkl_malloc(nrows * sizeof(double), 64);
  int*    nnzPerRow = (int*) mkl_malloc(nrows * sizeof(int), 64);
  int*    nonzeros = (int*) mkl_malloc(nrows*nnz_per_row * sizeof(int), 64);
  double* RsqContribs = (double*) mkl_malloc(nchunks * sizeof(double), 64);
  double* pApContribs = (double*) mkl_malloc(nchunks * sizeof(double), 64);
  double* Alpha = (double*) mkl_malloc(sizeof(double), 64);
  double* Beta = (double*) mkl_malloc(sizeof(double), 64);
  double* pAp = (double*) mkl_malloc(sizeof(double), 64);
  double* Rsq = (double*) mkl_malloc(sizeof(double), 64);
  double* RsqNextIter = (double*) mkl_malloc(sizeof(double), 64);
  double* b = (double*) mkl_malloc(nrows * sizeof(double), 64);
  for(int i = 0; i < nrows; i++){
    b[i] = i;
  }

  double** pChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);
  double** bChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);
  double** xChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);
  double** rChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);

  double** AChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);
  double** ApChunks = (double**) mkl_malloc(nchunks * sizeof(double*), 64);
  int**    nnzPerRowChunks = (int**) mkl_malloc(nchunks * sizeof(int*), 64);
  int**    nonzerosChunks = (int**) mkl_malloc(nchunks * sizeof(int*), 64);

  for (int i=0; i < nchunks; ++i){
    pChunks[i] = p + i*chunkSize;
    xChunks[i] = x + i*chunkSize;
    rChunks[i] = r + i*chunkSize;
    bChunks[i] = b + i*chunkSize;
    nnzPerRowChunks[i] = nnzPerRow + i*chunkSize;
    ApChunks[i] = Ap + i*chunkSize;
  }
  generate_problem_27pt(nx,ny,nz,chunkSize,A,nonzeros,nnzPerRow,AChunks,nonzerosChunks);

#pragma offload_transfer target(mic:0) \
  in(p : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(x : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(r : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(b : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(A : length(nrows*nnz_per_row) align(64) alloc_if(1) free_if(0)) \
  in(Ap : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(nnzPerRow : length(nrows) align(64) alloc_if(1) free_if(0)) \
  in(nonzeros : length(nrows*nnz_per_row) align(64) alloc_if(1) free_if(0)) \
  in(RsqContribs : length(nchunks) align(64) alloc_if(1) free_if(0)) \
  in(pApContribs : length(nchunks) align(64) alloc_if(1) free_if(0)) \
  in(Alpha : length(1) align(64) alloc_if(1) free_if(0)) \
  in(Beta : length(1) align(64) alloc_if(1) free_if(0)) \
  in(pAp : length(1) align(64) alloc_if(1) free_if(0)) \
  in(Rsq : length(1) align(64) alloc_if(1) free_if(0)) \
  in(RsqNextIter : length(1) align(64) alloc_if(1) free_if(0))

  Task* root = 0;
  root = initDag(cfg,
                 p, x, //full arrays
                 pApContribs, RsqContribs, //reduction arrays
                 Alpha, Beta, pAp, Rsq, RsqNextIter, //double pointers
                 AChunks, ApChunks, bChunks, pChunks, rChunks, xChunks, //double chunk arrays,
                 nnzPerRowChunks, nonzerosChunks //int chunk arrays
                );
  sch->run(root);

  return 0;
}

RegisterTest("cg", cg);
RegisterTest("cgprofiling", cgprofiling);

