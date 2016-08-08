#ifndef _lu_h_
#define _lu_h_

enum fxn_id {
  getrf_id,
  lu_gemm_id,
  trsm_left_id, 
  trsm_right_id, 
  dogetrf_id,
  dolu_gemm_id,
  dotrsm_left_id,
  dotrsm_right_id,
  getrfprofile_id,
  lu_gemmprofile_id,
  trsm_leftprofile_id,
  trsm_rightprofile_id,
};

void trsm_left(int k, int m, int size, int lda, double* L, double* A);
void trsm_right(int k, int m, int size, int lda, double* A, double* U);
void getrf(int k, int size, int lda, int* ipiv, double* A);
void lu_gemm(int m, int n, int k, int size, int lda, bool ltrans,  bool rtrans, double alpha, double beta, double* product,  double* left,  double* right);

#endif
