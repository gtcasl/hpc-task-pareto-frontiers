

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

