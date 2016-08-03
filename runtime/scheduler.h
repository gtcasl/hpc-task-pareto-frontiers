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
 private:
  virtual void runMaster(Task* root) = 0;

  static Scheduler* global;

  CPUList available_cores_;

  mic_device* mic_device_;

 protected:
  long long cumulative_power_;
  long num_power_samples_;
  uint32_t max_power_;
  double power_limit_; // I think the default of 350 is way above the TDP
  double available_power_; // amount of power that can be used at the moment

};

class SequentialScheduler : public Scheduler
{
  void runMaster(Task* root);
};

class BaselineScheduler : public Scheduler
{
  void runMaster(Task* root);
};

#endif

