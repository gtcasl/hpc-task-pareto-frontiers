#include <test.h>

typedef Buffer<double> DoublePtr;
typedef Buffer<double> IntPtr;

struct config
{
  int nchunks;
  int nnz_per_row;
  int nrows;
  int ncols;
};

enum fxn_id {
  ddot_id,
  start_id,
  daxpy_id,
  spmv_id
} test_ids;

void ddot(int chunk, int n, Buffer<double> a, Buffer<double> b, Buffer<double> result)
{
  double tmp 
  printf("Running ddot %p = %p * %p\n", c.buffer, a.buffer, b.buffer);
  for (int i=0; i < n; ++i){
    c[i] = a[i]*b[i];
  }
}

void start()
{
  //no op for organizational purposes
}

#define new_task(taskName, ...) \
  task(taskName, taskName_##id, std::make_tuple(__VA_ARGS__))

Task* initDag(config cfg, 
  DoublePtr A,
  DoublePtr Ax,
  DoublePtr p,
  DoublePtr r,
  DoublePtr x,
  const std::vector<DoublePtr>& pChunks,
  const std::vector<DoublePtr>& xChunks,
  const std::vector<DoublePtr>& rChunks,
  const std::vector<DoublePtr>& AChunks,
  const std::vector<DoublePtr>& AxChunks,
  const std::vector<DoublePtr>& nnzPerRowChunks,
{
  Task* root = task(start, start_id, std::make_tuple());
  std::vector<Task*> lastWavefront(cfg.nchunks);
  for (int i=0; i < cfg.nchunks; ++i){
    //have to start by computing Ax
    Task* compAx = new_task(spmv, cfg, AChunks[i], x, AxChunks[i]);
    Task* compR0 = new_task(daxpy, cfg, -1, bChunks[i], AxChunks[i]);
    Task* compP0 = new_task(memcpy, cfg, pChunks[i], rChunks[i]);
    compP0->dependsOn(compR0);
    compR0->dependsOn(compAx);
    lastWavefront[i] = compP0;
  }
  //now we enter the loop

  for (int i=0; i < cfg.nchunks; ++i){
    Task* compRsq   = new_task(ddot,  cfg, i, rsqChunks, rChunks[i], rChunks[i]);
    Task* compAp    = new_task(spmv,  cfg, AChunks[i], p, ApChunks[i]);
    Task* comp_pAp  = new_task(ddot,  cfg, i, pApChunks, ApChunks[i], pChunks[i]);
    Task* compAlpha = new_task(alpha, cfg, alpha, pApChunks, rsqChunks);
    Task* compX     = new_task(daxpy, cfg, alpha, xChunks[i], pChunks[i]);

    compAlpha->dependsOn(compRsq);
    comp_pAp->dependsOn(compAp);
    compX->dependsOn(compAlpha);
  }
};

int cg(int argc, char** argv)
{
  config cfg;
  cfg.nchunks = 2;
  cfg.nnz_per_row = 27;
  cfg.nrows = nx*ny*nz;
  cfg.ncols = cfg.nrows;

  int chunkSize = ncols / chunks;

  DoublePtr p(nrows);
  DoublePtr x(nrows);
  DoublePtr r(nrows);
  DoublePtr A(nrows*nnz_per_row);  
  IntPtr    nnzPerRow(nrows);
  IntPtr    nonzeros(nrows*nnz_per_row);
  IntPtr    rowOffsets(nrows);
  DoublePtr reduceRsq(nchunks);
  DoublePtr reduceBeta(nchunks);
  DoublePtr reduceAlpha(nchunks);
  DoublePtr reducepAp(nchunks);
  DoublePtr Beta(1);
  DoublePtr Alpha(1);
  DoublePtr RsqK(1);
  DoublePtr RsqKplus1(1);

  std::vector<DoublePtr> pChunks(nchunks);
  std::vector<DoublePtr> xChunks(nchunks);
  std::vector<DoublePtr> rChunks(nchunks);
  std::vector<DoublePtr> Achunks(nchunks);
  std::vector<IntPtr>    nnzPerRowChunks(nchunks);
  std::vector<IntPtr>    nonzerosChunks(nchunks);
  std::vector<IntPtr>    rowOffsets(nchunks);

  for (int i=0; i < nchunks; ++i){
    pChunks[0] = p.offset(i*chunkSize);
    xChunks[0] = x.offset(i*chunkSize);
    rChunks[0] = r.offset(i*chunkSize);
    AChunks[0];
  }

  //set up the problem



  return 0;
}

