#include "index.h"
#include <cstdio>
#include <iostream>

void
generate_problem_27pt(
  int nx, int ny, int nz,
  int* nrows_ptr,
  int** nnzPerRow_ptr,
  int*** nonzerosInRow_ptr,
  double*** matrixRows_ptr,
  int** indices_ptr,
  double** matrix_ptr)
{
  int nnz_per_row = 27;
  int nrows = nx*ny*nz;
  int ncols = nrows;
  double* vector = new double[nrows];
  int* nnzPerRow = new int[nrows];
  int* nonzeros = new int[nrows*nnz_per_row];
  double* matrix = new double[nrows*nnz_per_row];

  double* matrixValPtr = matrix;
  int* indexValPtr = nonzeros;
  
  double** matrixRows = new double*[nrows];
  int** nonzerosInRow = new int*[nrows];

  for (int rx=0; rx < nx; rx++){
    int cxStart = max(rx-1,0);
    int cxStop = min(rx+1,nx-1);
    for (int ry=0; ry < ny; ry++){
      int cyStart = max(ry-1,0);
      int cyStop = min(ry+1,ny-1);
      for (int rz=0; rz < nz; ++rz){
        int czStart = max(rz-1,0);
        int czStop = min(rz+1,nz-1);
        int row = index(rx,ry,rz,nx,ny,nz);
        int nnzInRow = 0;
        nonzerosInRow[row] = indexValPtr;
        matrixRows[row] = matrixValPtr;
        for (int cx=cxStart; cx <= cxStop; ++cx){
          for (int cy=cyStart; cy <= cyStop; ++cy){
            for (int cz=czStart; cz <= czStop; ++cz, ++matrixValPtr, ++indexValPtr, ++nnzInRow){
              int col = index(cx,cy,cz,nx,ny,nz);
              *indexValPtr = col;
              if (row==col){
                *matrixValPtr = 26.0;
              } else {
                *matrixValPtr = -1.0;
              }
            }
          }
        }
        nnzPerRow[row] = nnzInRow;
      }
    }
  }

  *nnzPerRow_ptr = nnzPerRow;
  *nonzerosInRow_ptr = nonzerosInRow;
  *matrixRows_ptr = matrixRows;
  *indices_ptr = nonzeros;
  *matrix_ptr = matrix;
  *nrows_ptr = nrows;
}

