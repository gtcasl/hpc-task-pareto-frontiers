#ifndef _test_runtime_scheduler_h_
#define _test_runtime_scheduler_h_

#define heisenbug printf("heisenbug %s:%d\n", __FILE__, __LINE__)

#include <vector>
#include <task.h>
#include <unordered_set>
#ifndef no_miclib
#include <miclib.h>
#endif
#include <signal.h>

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

  void resetIter();

  void*
  relocatePointer(size_t offset){
    void* ret = ((char*)mmap_buffer_) + offset;
    return ret;
  }

  int numAvailableCores() const {
    return available_cores_.size();
  }

  void returnCpu(int cpu);

  int claimCpu();

  virtual ~Scheduler(){
    global = 0;
#ifndef no_miclib
    mic_close_device(mic_device_);
#endif
  }

  uint32_t readMICPoweruW() const;
  static void overflow(int signum, siginfo_t*, void*);
  double getTime() const;

 protected:
  Scheduler() :
    mmap_buffer_(0),
    mmap_size_(0),
    total_buffer_size_(0),
    next_copy_(0),
    cumulative_power_(0),
    num_power_samples_(0),
    max_power_(0),
    power_limit_(1000),
    available_power_(power_limit_)
  {
    if (global){
      fprintf(stderr, "only allowed one instance of a scheduler at a time\n");
      abort();
    }
    global = this;
    // initialize power measurement
#ifndef no_miclib
    if(mic_open_device(&mic_device_, 0) != E_MIC_SUCCESS){
      std::cerr << "Error: Unable to open MIC" << std::endl;
      abort();
    }
#endif
  }
  long long cumulative_power_;
  long num_power_samples_;
  uint32_t max_power_;
  uint32_t power_limit_; // I think the default of 350 is way above the TDP
  uint32_t available_power_; // amount of power that can be used at the moment

 protected:
  void runSerial(Task* task);

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

  std::unordered_set<int> available_cores_;

#ifndef no_miclib
  struct mic_device* mic_device_;
#endif
};

template <class T>
class Buffer : public BufferBase {
 public:
  Buffer(int n) {
    nelems = n;
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
    Buffer<T> newbuf(nelems, mmap_offset + off*sizeof(T));
    T* relocatedBuf = (T*) Scheduler::global->relocatePointer(mmap_offset);
    newbuf.buffer = (relocatedBuf + off);
    return newbuf;
  }

  const T&
  operator[](int idx) const {
    const T* buf = (const T*) buffer;
    return buf[idx];
  }

};

template <class T>
std::ostream& 
operator<<(std::ostream& os, const Buffer<T>& buf){
  os << "Buffer: offset=" << buf.mmap_offset << " ptr=" << (void*) buf.buffer << " nelems=" << buf.nelems;
  return os;
}

class BasicScheduler : public Scheduler
{
 public:
  void runMaster(Task* root);
};

class AdvancedScheduler : public Scheduler
{
 public:
  void runMaster(Task* root);
};



#endif

