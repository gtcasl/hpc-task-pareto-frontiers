#ifndef _test_runtime_scheduler_h_
#define _test_runtime_scheduler_h_

#include <vector>
#include <task.h>
#include <unordered_set>
#include <signal.h>

struct mic_device;

double getTime();

struct CPUList{
  int num_available;
  std::vector<std::pair<int,bool>> cpus; // pair(cpu_id, available?)
  CPUList() {}
  CPUList(int na) : num_available{na} {
    for(int i = 0; i < num_available; i++){
      auto cpuid = i % 57 * 4 + i / 57;
      cpus.emplace_back(cpuid, true);
    }
  }
};

class Scheduler
{

 public:
  void run(Task* root);

  void finalize();

  virtual void init(int argc, char** argv);

  int numAvailableCores() const {
    return available_cores_.num_available;
  }

  void returnCpu(int cpu);

  int claimCpu();

  virtual ~Scheduler();

  uint32_t readMICPoweruW() const;
  static void overflow(int signum, siginfo_t*, void*);

 protected:
  Scheduler();
  virtual void tick() = 0;
 private:

  static Scheduler* global;

  CPUList available_cores_;

  mic_device* mic_device_;

 protected:
  long long cumulative_power_;
  long num_power_samples_;
  uint32_t max_power_;
  double power_limit_; // I think the default of 350 is way above the TDP
  double current_power_; // estimated total of current power usage
  bool do_profiling_;
  std::ofstream logfile_;
  std::list<pthread_t> runningTasks_;
  std::list<Task*> pendingTasks_;
  std::map<Task*, std::unordered_set<int> > taskCpuAssignments_;
  int tick_number_;
  double estimated_max_power_;

  void launchTask(Task* task, int nthreads);

  void log(){
    logfile_ << std::endl;
  }

  template<class T, class U, class... Ts>
  void log(const T& t, const U& u, const Ts&... args){
    logfile_ << t << "=" << u << ",";
    log(args...);
  }
};

class SequentialScheduler : public Scheduler
{
  public:
    SequentialScheduler(bool do_prof){
      do_profiling_ = do_prof;
    }
  protected:
    void tick();
};

class FairScheduler : public Scheduler
{
  protected:
    void tick();
};

class CoreConstrainedScheduler : public Scheduler
{
  protected:
    void tick();
};

class PowerAwareScheduler : public Scheduler
{
  protected:
    void tick();
    virtual void leverageSlack(std::map<Task*, std::vector<ParetoPoint>::iterator>&) {}
};

class SlackAwareScheduler : public PowerAwareScheduler
{
  protected:
    void leverageSlack(std::map<Task*, std::vector<ParetoPoint>::iterator>& task_assignments);
};

#endif

