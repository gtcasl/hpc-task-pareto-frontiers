#ifndef _test_runtime_scheduler_h_
#define _test_runtime_scheduler_h_

#include <vector>
#include <task.h>

class BufferBase
{
 public:
  size_t mmap_offset;
  int nelems;
  void* buffer;
};

class Scheduler
{
  template <class T> friend class Buffer;

 public:
  void run(Task* root);

  void finalize();

  void stop();

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

  void nextIter();

  void*
  relocatePointer(size_t offset){
    void* ret = ((char*)mmap_buffer_) + offset;
    return ret;
  }

  virtual ~Scheduler(){
    global = 0;
  }

 protected:
  Scheduler() :
    mmap_buffer_(0),
    mmap_size_(0),
    total_buffer_size_(0),
    next_copy_(0) {
    if (global){
      fprintf(stderr, "only allowed one instance of a scheduler at a time\n");
      abort();
    }
    global = this;
  }

 private:
  void addNeededBuffer(BufferBase* buf, size_t size);

  /**
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
  */

  void assignBuffer(BufferBase* buf);
  /**
  template<class... Ts>
  void assignBuffers(Buffer<Ts>&... bufs){
    // see the note for addNeededBuffers for explanation of variadic template
    int dummy[] = { 0, (assignBuffer(bufs), 0)... };
  }
  */

  double getTime() const;


  void runWorker();
  void terminateWorkers();
  virtual void runMaster(Task* root) = 0;

  void* mmap_buffer_;
  size_t mmap_size_;
  size_t total_buffer_size_;
  int next_copy_;
  int ncopies_;
  std::list<BufferBase*> buffers_;

  static Scheduler* global;

  int nworkers_;
  int rank_;
  int nproc_;
};

template <class T>
class Buffer : public BufferBase {
 public:
  Buffer(int n) {
    nelems = 0;
    mmap_offset = 0;
    buffer = 0;
    Scheduler::global->addNeededBuffer(this, n*sizeof(T));
  }

  Buffer(int n, int offset)
  {
    mmap_offset = offset;
    buffer = 0;
    nelems = n;
  }

  Buffer(const Buffer<T>& t)
  {
    mmap_offset = t.mmap_offset;
    nelems = t.nelems;
    buffer = (T*) Scheduler::global->relocatePointer(mmap_offset);
  }

  operator T*() {
    return (T*) buffer;
  }

  T&
  operator[](int idx){
    T* buf = (T*) buffer;
    return buf[idx];
  }

  T&
  operator*(){
    T* buf = (T*) buffer;
    return *buf;
  }

  Buffer<T>
  offset(int off) const {
    return Buffer<T> (nelems, mmap_offset + off*sizeof(T));
  }

  const T&
  operator[](int idx) const {
    const T* buf = (const T*) buffer;
    return buf[idx];
  }

};

class BasicScheduler : public Scheduler
{
 public:
  void runMaster(Task* root);
};



#endif

