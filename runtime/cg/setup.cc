#include <iostream>

inline int
index(int ix, int iy, int iz, int, int ny, int nz){
  return ix*ny*nz + iy*nz + iz;
}

inline int
max(int l, int r){
  return l < r ? r : l;
}

inline int
min(int l, int r){
  return l < r ? l : r;
}

void
generate_problem_27pt(
  int nx, int ny, int nz,
  int chunkSize,
  double* A,
  int* nonzeros,
  int* nnzPerRow,
  double** AChunks,
  int** nonzerosChunks
)
{

  int offset = 0;
  int nextChunk = 0;
//#pragma omp parallel for
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
        if (row % chunkSize == 0){
          AChunks[nextChunk] = A + offset;
          nonzerosChunks[nextChunk] = nonzeros + offset;
          ++nextChunk;
        }
        for (int cx=cxStart; cx <= cxStop; ++cx){
          for (int cy=cyStart; cy <= cyStop; ++cy){
            for (int cz=czStart; cz <= czStop; ++cz, ++offset, ++nnzInRow){
              int col = index(cx,cy,cz,nx,ny,nz);
              nonzeros[offset] = col;
              if (row==col){
                A[offset] = 26.0;
              } else {
                A[offset] = -1.0;
              }
            }
          }
        }
        if (nnzInRow > 27) abort();
        nnzPerRow[row] = nnzInRow;
      }
    }
  }
}

