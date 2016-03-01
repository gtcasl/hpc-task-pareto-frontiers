#include <sched.h>
#include <type_traits>
#include <mpi.h>
#include <tuple>
#include <map>
#include <list>
#include <set>
#include <utility>

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

  int typeID() const { return typeID_; }

  void run(int worker);
  
  void setup();
  
  void notifyDone();

  int worker() const {
    return worker_;
  }

  void dependsOn(Task* t, int bytes = 0){
    preconditions_[t] += bytes;
    t->listeners_.insert(this);
  }

  int clearDependency(Task* t){
    preconditions_.erase(t);
    return preconditions_.size();
  }

  void clearListeners(std::list<Task*>& ready);
  
 protected:
  Task(int mySize, int typeID);

 private:
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

};

class TaskRunner {
 public:
  virtual void run(Task* t) = 0;

  static void store(int id, TaskRunner* r){
    runners_[id] = r;
  }

  static TaskRunner* get(int id){
    return runners_[id];
  }

 private:
  static std::map<int, TaskRunner*> runners_;
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


template <typename Fxn, typename First, typename ...Args>
class TaskRunner_impl : public TaskRunner
{
  typedef Task_tmpl<Fxn,First,Args...> TaskType;
 public:
  void run(Task* t)
  {
    auto theTask = static_cast<TaskType*>(t);
    theTask->setup();
    theTask->run();
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
registerFunction(typename myFxnTraits::Fxn f, int id){
  impl::FunctionMap<typename myFxnTraits::Fxn>::store(id, f);
  typedef typename myFxnTraits::Runner myRunner;
  TaskRunner::store(id, new myRunner);
  return 0;
}

#define RegisterTask(fxn,id,ret,...) \
  typedef impl::FxnTraits<ret,__VA_ARGS__> type_traits_##fxn; \
  int ignore_##fxn = registerFunction<type_traits_##fxn>(fxn,id);
  

template <typename Fxn, typename... Args>
Task*
task(Fxn f, int id, const std::tuple<Args...>& t){
  return new impl::Task_tmpl<Fxn,Args...>(t,id);
}

