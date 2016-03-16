#ifndef _test_runtime_scheduler_h_
#define _test_runtime_scheduler_h_

#include <vector>
#include "task.h"

template <class T>
class Buffer {
 public:
  size_t mmap_offset;
  int nelems;
  T* buffer;

  Buffer(int n) : mmap_offset(0), buffer(0), nelems(n) {}

  Buffer(const Buffer<T>& t);

  Buffer(Buffer<T>&&) = delete;
  
  T&
  operator[](int idx){
    return buffer[idx];
  }

  const T&
  operator[](int idx) const {
    return buffer[idx];
  }

};

class Scheduler
{
 public:
  void run(Task* root);
  
  virtual void init(int argc, char** argv);

  static const int terminate_tag = 42;

  int rank() const { 
    return rank_;
  }

  int nworkers() const {
    return nworkers_;
  }

  int nproc() const {
    return nproc_;
  }

  void allocateHeap(int ncopies);

  void deallocateHeap();

  template <class T>
  void addNeededBuffer(Buffer<T>& buf){
    size_t size = buf.nelems * sizeof(T);
    buf.mmap_offset = total_buffer_size_;
    if (size % 4096){
      size_t extra = 4096 - size % 4096;
      size += extra;
    }
    total_buffer_size_ += size;
  }

  template <class... Ts>
  void addNeededBuffers(Buffer<Ts>&... bufs){
    // To call a function on each element of the parameter pack, 
    // we need to play some tricks. First, we make a dummy array and
    // take advantage of initializer list expansion for parameter packs.
    // Then, since the function we call returns void, we have to force it
    // to return a number, using the comma operator. Finally, since this 
    // function can be called with zero parameters, we need to make sure 
    // that the dummy array is at least one element long.
    int dummy[] = { 0, (addNeededBuffer(bufs), 0)... };
  }

  template <class T>
  void assignBuffer(Buffer<T>& buf){
    if (next_copy_ % ncopies_ == 0){
      buf.mmap_offset = buf.mmap_offset % total_buffer_size_;
    }
    else {
      buf.mmap_offset += total_buffer_size_;
    }
  }

  template<class... Ts>
  void assignBuffers(Buffer<Ts>&... bufs){
    // see the note for addNeededBuffers for explanation of variadic template
    int dummy[] = { 0, (assignBuffer(bufs), 0)... };
  }

  void nextIter() {
    next_copy_++;
  }

  static void*
  relocatePointer(size_t offset){
    void* ret = ((char*)mmap_buffer_) + offset;
    return ret;
  }

 protected:
  Scheduler() : total_buffer_size_(0), next_copy_(0) {}
 
 private:
  void runWorker();
  void terminateWorkers();
  virtual void runMaster(Task* root) = 0;

  static void* mmap_buffer_;
  static size_t mmap_size_;
  int nworkers_;
  int rank_;
  int nproc_;
  int next_copy_;
  int ncopies_;

  size_t total_buffer_size_;

};


class BasicScheduler : public Scheduler
{
 public:
  void runMaster(Task* root);
};

template <class T>
Buffer<T>::Buffer(const Buffer<T>& t){
  mmap_offset = t.mmap_offset;
  nelems = t.nelems;
  buffer = (T*) Scheduler::relocatePointer(mmap_offset);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //printf("relocated pointer at offset %lu to %p for n=%d : %p->%p on rank %d\n", 
  //  mmap_offset, buffer, nelems, &t, this, rank);
}

#endif

