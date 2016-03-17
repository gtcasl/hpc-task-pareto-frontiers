#include "index.h"
#include <cstdio>
#include <iostream>

void
generate_problem_27pt(
  int nrows,
  int nnz_per_row,
  int* nnzPerRow,
  int* rowOffsets,
  double* matrix,
  int* nonzeroes)
{
  int nrows = nx*ny*nz;
  int ncols = nrows;

  double* matrixValPtr = matrix;
  int* indexValPtr = nonzeros;
  
  int offset = 0;
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
        rowOffsets[row] = offset;
        for (int cx=cxStart; cx <= cxStop; ++cx){
          for (int cy=cyStart; cy <= cyStop; ++cy){
            for (int cz=czStart; cz <= czStop; ++cz, ++offset, ++indexValPtr, ++nnzInRow){
              int col = index(cx,cy,cz,nx,ny,nz);
              *indexValPtr = col;
              if (row==col){
                matrix[offset] = 26.0;
              } else {
                matrix[offset] = -1.0;
              }
            }
          }
        }
        nnzPerRow[row] = nnzInRow;
      }
    }
  }
}

