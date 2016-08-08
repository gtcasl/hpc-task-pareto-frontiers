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

#endif
