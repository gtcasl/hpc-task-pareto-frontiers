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
#include <cassert>
#include <tuple>

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

#define new_mutating_task(taskName, ...) \
  make_task(taskName, #taskName, taskName##_id, std::make_tuple(__VA_ARGS__), true)

#define new_task(taskName, ...) \
  make_task(taskName, #taskName, taskName##_id, std::make_tuple(__VA_ARGS__), false)


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

  double estimateMakespan(int nthreads) const;

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

  void runProfiling();

  bool canProfile() const {
    return !isMutating;
  }

  int getIters() const {
    return numIters;
  }

 protected:
  Task(int typeID, bool isMut);
  int numIters;

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
  bool isMutating;
};

struct ParetoPoint{
  int nthreads;
  double time;
  double power;
  ParetoPoint(int n, double t, double p) : nthreads{n}, time{t}, power{p} {}
};

extern std::map<int, const char*> Names;
extern std::map<int, std::vector<double> > Times;
extern std::map<int, std::vector<double> > Powers;
extern std::map<int, std::vector<ParetoPoint>> Paretos;

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
  Task_tmpl(Fxn f, const Tuple& params, int typeID, bool isMut) :
    Task(typeID, isMut),
    params_(params), func_(f)
  {}

  void do_run() {
    numIters++;
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

std::vector<ParetoPoint>
make_pareto(const std::vector<double>& time,
            const std::vector<double>& power);

template <typename Fxn, typename... Args>
Task*
make_task(Fxn f, const char* name, int id, const std::tuple<Args...>& t, bool isMut){
  if(Names.find(id) == std::end(Names)){
    // insert all the stuff into the static maps
    Names[id] = name;
    std::vector<double> times(NUM_THREADS + 1, std::numeric_limits<double>::infinity());
    std::vector<double> powers(NUM_THREADS + 1, std::numeric_limits<double>::infinity());
    powers[0] = 0.0; // default case is infinite power, but no threads means no power
    std::string sname = name;
    sname += ".csv";
    std::ifstream ifs{sname};
    if(ifs.good()){
      std::string line;
      std::getline(ifs, line); // kill header
      while(std::getline(ifs, line)){
        std::stringstream stst{line};
        std::string elem;
        std::getline(stst, elem, ','); // nthreads
        int nthreads = std::stoi(elem);
        assert(nthreads < NUM_THREADS && "Error: profile uses more than the available number of threads");
        std::getline(stst, elem, ','); // energy
        std::getline(stst, elem, ','); // time
        times[nthreads] = std::stod(elem);
        std::getline(stst, elem, ','); // power
        // subtract out the baseline for keeping the system running
        powers[nthreads] = std::stod(elem) - IDLE_POWER;
        assert(powers[nthreads] >= 0.0 && "Error: power is lower than the idle power");
        std::getline(stst, elem, ','); // speedup
      }
      Paretos[id] = make_pareto(times, powers);
    } else {
      // make it so only a single thread is assigned whenever a csv doesn't exist
      std::fill(std::begin(times)+1, std::end(times), 0);;
      std::fill(std::begin(powers)+1, std::end(powers), 0);;
      Paretos[id].emplace_back(1,0,0);
      //times[1] = 0;
      //powers[1] = 0;
    }
    Times[id] = times;
    Powers[id] = powers;
    /*
    std::cout << "Pareto for " << name << "\n";
    for(const auto& e : Paretos[id]){
      std::cout << e.nthreads << "\t" << e.time << "\t" << e.power << "\n";
    }
    */
  }
  return new impl::Task_tmpl<Fxn,Args...>(f,t,id,isMut);
}

#endif
