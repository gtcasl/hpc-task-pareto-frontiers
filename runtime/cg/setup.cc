#include <cg/types.h>

typedef Buffer<double> DoublePtr;
typedef Buffer<double> DoubleArray;
typedef Buffer<int> IntArray;
typedef Buffer<IntArray> IntChunkArray;
typedef Buffer<DoubleArray> DoubleChunkArray;

inline int
index(int ix, int iy, int iz, int nx, int ny, int nz){
  return iz*nx*ny + iy*nx + ix;
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
  DoubleArray A,
  IntArray nonzeros,
  IntArray nnzPerRow,
  DoubleChunkArray bChunks,
  DoubleChunkArray xChunks,
  DoubleChunkArray AChunks,
  IntChunkArray    nonzerosChunks
)
{
  int nrows = nx*ny*nz;
  int ncols = nrows;


  int offset = 0;
  int nextChunk = 0;
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
          AChunks[nextChunk] = A.offset(offset);
          nonzerosChunks[nextChunk] = nonzeros.offset(offset);
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
        nnzPerRow[row] = nnzInRow;
      }
    }
  }
}

