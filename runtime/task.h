#ifndef _test_runtime_task_h_
#define _test_runtime_task_h_

#include <iostream>
#include <sched.h>
#include <type_traits>
#include <mpi.h>
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
  make_task(taskName, taskName##_id, std::make_tuple(__VA_ARGS__))


#ifdef DEP_DEBUG
#define dep_debug(...) printf(__VA_ARGS__)
#else
#define dep_debug(...) 
#endif

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

  bool checkDone();

  void waitDone();

  cpu_set_t*
  cpuMask() {
    return &cpumask_;
  }

  int
  cpuCount() const {
    return CPU_COUNT(&cpumask_);
  }

  int typeID() const { return typeID_; }

  void run(int worker, int start_tick);
  
  void setup();
  
  void notifyDone();

  int worker() const {
    return worker_;
  }

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

  double getStartTime() const {
    return start_;
  }

  int getStartTick() const { return start_tick_; }

  int getNumThreads() const;

  double estimateTime() const;

  int getNumListeners() const { return listeners_.size(); }
  
 protected:
  Task(int mySize, int typeID);

 private:
  double getTime() const;

  static const int enqueue_tag = 100;
  static const int done_tag = 101;
  int mySize_;
  int cpu_state_;
  int nthread_;
  int typeID_;
  int worker_;
  bool done_;
  cpu_set_t cpumask_;
  MPI_Request done_request_;
  int rc_;
  /**
   * Key: Task ID
   * Value: Number of bytes in the dependency
   */
  std::map<Task*,int> preconditions_;
  /**
   * The tasks that require this task to complete before running
   */
  std::set<Task*> listeners_;
  double start_;
  int start_tick_;
};

class TaskRunner {
 public:
  virtual void disableProfiling(){}
  virtual void run(Task* t, int size) = 0;

  static void store(int id, TaskRunner* r, const char* name){
    runners_[id] = r;
    names_[id] = name;
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
    times_[id] = times;
    powers_[id] = powers;
    auto min = std::min_element(std::begin(times), std::end(times));
    min_times_[id] = *min;
    min_threads_[id] = std::distance(std::begin(times), min);
  }

  static TaskRunner* get(int id){
    return runners_[id];
  }

  static const char* get_name(int id){
    return names_[id];
  }

  static const std::vector<double>& get_times(int id){
    return times_[id];
  }

  static const std::vector<double>& get_powers(int id){
    return powers_[id];
  }

  static double get_min_time(int id){
    return min_times_[id];
  }

  static int get_min_threads(int id){
    return min_threads_[id];
  }

  static int get_next_least_powerful_num_threads(int id, int cur_num_threads);

 private:
  static std::map<int, TaskRunner*> runners_;
  static std::map<int, const char*> names_;
  static std::map<int, std::vector<double> > times_;
  static std::map<int, double> min_times_;
  static std::map<int, int> min_threads_;
  static std::map<int, std::vector<double> > powers_;
};

namespace impl {

template <typename Fxn>
class FunctionMap
{
 public:
  static Fxn get(int id){
    return functions_[id];
  }
  static void store(int id, Fxn f){
    functions_[id] = f;
  }
 private:
  static std::map<int,Fxn> functions_;
};

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
  Task_tmpl(const Tuple& params, int typeID) :
    Task(sizeof(Task_tmpl), typeID),
    params_(params)
  {
  }

  void run() {
    callFunc(typename gens<sizeof...(Args)>::type());
  }

 private:
  Tuple params_;

  template<int ...S>
  void callFunc(seq<S...>) {
    Fxn f = FunctionMap<Fxn>::get(typeID());
    return f(std::get<S>(params_) ...);
  }
};

extern bool done;

void wakeup(int signum, siginfo_t*, void*);

template <typename Fxn, typename First, typename ...Args>
class TaskRunner_impl : public TaskRunner
{
  typedef Task_tmpl<Fxn,First,Args...> TaskType;
 public:
  void disableProfiling(){doProfiling = false;}
  bool doProfiling;
  TaskRunner_impl(bool profile) : doProfiling{profile}{}
  void run(Task* t, int size){
    if(doProfiling){
      runProfiling(t, size);
    } else {
      runNoProfiling(t, size);
    }
  }
  void runNoProfiling(Task* t, int size)
  {
    if (sizeof(TaskType) != size){
      std::cerr << "type mismatch between task received and runner on ID " 
        << t->typeID() << std::endl;
      abort();
    }
    auto theTask = static_cast<TaskType*>(t);

    theTask->setup();
    theTask->run();
    theTask->notifyDone();
  }

  void runProfiling(Task* t, int size)
  {
    if (sizeof(TaskType) != size){
      std::cerr << "type mismatch between task received and runner on ID " 
        << t->typeID() << std::endl;
      abort();
    }
    auto theTask = static_cast<TaskType*>(t);

    done = false;
    struct sigaction sa;
    sa.sa_sigaction = wakeup;
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

    theTask->setup();
    while(!done){
      theTask->run();
    }
    done = false;
    theTask->notifyDone();
  }
};

template <typename R, typename... Args>
class FxnTraits
{
 public:
  typedef R (*Fxn)(Args...);
  typedef impl::TaskRunner_impl<Fxn,Args...> Runner;
  typedef impl::Task_tmpl<Fxn,Args...> Task;
};

}


template <typename myFxnTraits>
int
registerFunction(typename myFxnTraits::Fxn f, int id, const char* name, bool profile){
  impl::FunctionMap<typename myFxnTraits::Fxn>::store(id, f);
  typedef typename myFxnTraits::Runner myRunner;
  TaskRunner::store(id, new myRunner(profile), name);
  return 0;
}

#define RegisterMutatingTask(fxn,ret,...) \
  typedef impl::FxnTraits<ret,__VA_ARGS__> type_traits_##fxn; \
  int ignore_##fxn = registerFunction<type_traits_##fxn>(fxn,fxn##_id,#fxn,false); \
  (void)ignore_##fxn;
  
#define RegisterTask(fxn,ret,...) \
  typedef impl::FxnTraits<ret,__VA_ARGS__> type_traits_##fxn; \
  int ignore_##fxn = registerFunction<type_traits_##fxn>(fxn,fxn##_id,#fxn,true); \
  (void)ignore_##fxn;
  

template <typename Fxn, typename... Args>
Task*
make_task(Fxn, int id, const std::tuple<Args...>& t){
  return new impl::Task_tmpl<Fxn,Args...>(t,id);
}

template <class T> std::map<int,T> impl::FunctionMap<T>::functions_;

#endif
