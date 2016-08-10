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
  virtual void runMaster(Task* root) = 0;
 private:

  static Scheduler* global;

  CPUList available_cores_;

  mic_device* mic_device_;

 protected:
  long long cumulative_power_;
  long num_power_samples_;
  uint32_t max_power_;
  double power_limit_; // I think the default of 350 is way above the TDP
  double available_power_; // amount of power that can be used at the moment
  bool do_profiling_;

  std::ofstream logfile_;
  void log(){
    logfile_ << std::endl;
  }

  template<class T, class U, class... Ts>
  void log(const T& t, const U& u, const Ts&... args){
    logfile_ << t << "=" << u << ",";
    log(args...);
  }
};

class FairScheduler : public Scheduler
{
  protected:
    void runMaster(Task* root);
};

class SequentialScheduler : public Scheduler
{
  protected:
    void runMaster(Task* root);
};

class SimpleScheduler : public Scheduler
{
  protected:
    void runMaster(Task* root);
};

class AdvancedScheduler : public Scheduler
{
  protected:
    void runMaster(Task* root);
};

class NonParetoScheduler : public Scheduler
{
  protected:
    void runMaster(Task* root);
};

class ProfilingScheduler : public SequentialScheduler
{
  protected:
    void runMaster(Task* root);
};

class BaselineScheduler : public AdvancedScheduler
{
  protected:
    void runMaster(Task* root);
};

#endif

