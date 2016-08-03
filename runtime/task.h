#ifndef _test_runtime_task_h_
#define _test_runtime_task_h_

#include <iostream>
#include <sched.h>
#include <type_traits>
#include <tuple>
#include <map>
#include <list>
#include <set>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <signal.h>

#ifdef no_affinity
struct cpu_set_t {};
#define sched_setaffinity(...)
#define CPU_ZERO(...)
#define CPU_SET(...)
#define CPU_COUNT(x)  1
#endif

#define NUM_THREADS 228
// This number was found empirically using the linear regression model from
// the energy scaling work.
#define IDLE_POWER 106.9

#define new_task(taskName, ...) \
  make_task(taskName, #taskName, taskName##_id, std::make_tuple(__VA_ARGS__))


#ifdef DEP_DEBUG
#define dep_debug(...) printf(__VA_ARGS__)
#else
#define dep_debug(...) 
#endif

namespace impl{
extern bool done;
void wakeup(int signum, siginfo_t*, void*);
}

class Task {
 public:
  void setCPUState(int s){
    cpu_state_ = s;
  }

  int CPUState() const {
    return cpu_state_; 
  }

  void setNumThreads(int n) {
    nthread_ = n;
  }

  void addCpu(int c);

  void addCpus(int cOffset, int numCores);

  cpu_set_t*
  cpuMask() {
    return &cpumask_;
  }

  int
  cpuCount() const {
    return CPU_COUNT(&cpumask_);
  }

  int typeID() const { return typeID_; }

  virtual void do_run() = 0;
  
  void setup();
  
  // you are allowed to pass a null task for 
  // algorithmic simplicity - it just gets ignored
  void dependsOn(Task* t, int bytes = 0){
    if (t){
      preconditions_[t] += bytes;
      t->listeners_.insert(this);
      dep_debug("Adding dependency %p for task %p - join counter = %d\n",
        t, this, preconditions_.size());
    }
  }

  int clearDependency(Task* t){
    preconditions_.erase(t);
    dep_debug("Clearing dependency %p for task %p - join counter = %d\n",
      t, this, preconditions_.size());
    return preconditions_.size();
  }

  void clearListeners(std::list<Task*>& ready);

  int getNumThreads() const;

  double estimateTime() const;

  int getNumListeners() const { return listeners_.size(); }

  void enableProfiling(){doProfiling = true;}

  void run(){
    if(doProfiling){
      runProfiling();
    } else {
      runNoProfiling();
    }
  }

  void runNoProfiling()
  {
    setup();
    do_run();
  }

  void runProfiling()
  {
    impl::done = false;
    struct sigaction sa;
    sa.sa_sigaction = impl::wakeup;
    sa.sa_flags = SA_SIGINFO;
    if(sigaction(SIGALRM, &sa, nullptr) != 0){
      std::cerr << "Error: unable to set up signal handler" << std::endl;
      abort();
    }
    struct itimerval work_time;
    work_time.it_value.tv_sec = 0;
    work_time.it_value.tv_usec = 500000; // sleep for 500 ms
    work_time.it_interval.tv_sec = 0;
    work_time.it_interval.tv_usec = 0;
    setitimer(ITIMER_REAL, &work_time, NULL);

    setup();
    while(!impl::done){
      do_run();
    }
    impl::done = false;
  }

 protected:
  Task(int typeID);

 private:
  int cpu_state_;
  int nthread_;
  int typeID_;
  cpu_set_t cpumask_;
  cpu_set_t full_mask_;
  /**
   * Key: Task ID
   * Value: Number of bytes in the dependency
   */
  std::map<Task*,int> preconditions_;
  /**
   * The tasks that require this task to complete before running
   */
  std::set<Task*> listeners_;
  bool doProfiling;
};

extern std::map<int, const char*> Names;
extern std::map<int, std::vector<double> > Times;
extern std::map<int, double> MinTimes;
extern std::map<int, int> MinThreads;
extern std::map<int, std::vector<double> > Powers;

int get_next_least_powerful_num_threads(int id, int cur_num_threads);

namespace impl {

template<int ...>
struct seq { };

template<int N, int ...S>
struct gens : gens<N-1, N-1, S...> { };

template<int ...S>
struct gens<0, S...> {
  typedef seq<S...> type;
};

template <typename Fxn, typename ...Args>
class Task_tmpl : public Task
{
  typedef std::tuple<Args...> Tuple;
 public:
  Task_tmpl(Fxn f, const Tuple& params, int typeID) :
    Task(typeID),
    params_(params), func_(f)
  {}

  void do_run() {
    callFunc(typename gens<sizeof...(Args)>::type());
  }

 private:
  Tuple params_;
  Fxn func_;

  template<int ...S>
  void callFunc(seq<S...>) {
    return func_(std::get<S>(params_) ...);
  }
};
}

template <typename Fxn, typename... Args>
Task*
make_task(Fxn f, const char* name, int id, const std::tuple<Args...>& t){
  if(Names.find(id) == std::end(Names)){
    // insert all the stuff into the static maps
    Names[id] = name;
    std::vector<double> times(NUM_THREADS + 1, 0.0);
    times[0] = std::numeric_limits<double>::infinity();
    std::vector<double> powers(NUM_THREADS + 1, 0.0);
    std::string sname = name;
    sname += ".csv";
    std::ifstream ifs{sname};
    if(ifs.good()){
      std::string line;
      std::getline(ifs, line); // kill header
      for(int i = 1; i < NUM_THREADS + 1; ++i){
        std::getline(ifs, line); 
        std::stringstream stst{line};
        std::string elem;
        std::getline(stst, elem, ','); // nthreads
        std::getline(stst, elem, ','); // energy
        std::getline(stst, elem, ','); // time
        times[i] = std::stod(elem);
        std::getline(stst, elem, ','); // power
        // subtract out the baseline for keeping the system running
        powers[i] = std::stod(elem) - IDLE_POWER;
        std::getline(stst, elem, ','); // speedup
      }
    }
    Times[id] = times;
    Powers[id] = powers;
    auto min = std::min_element(std::begin(times), std::end(times));
    MinTimes[id] = *min;
    MinThreads[id] = std::distance(std::begin(times), min);
  }
  return new impl::Task_tmpl<Fxn,Args...>(f,t,id);
}

#endif
