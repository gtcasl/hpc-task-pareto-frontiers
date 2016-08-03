#ifndef _test_runtime_scheduler_h_
#define _test_runtime_scheduler_h_

#include <vector>
#include <task.h>
#include <unordered_set>
#include <miclib.h>
#include <signal.h>

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

  virtual ~Scheduler(){
    global = 0;
    mic_close_device(mic_device_);
  }

  uint32_t readMICPoweruW() const;
  static void overflow(int signum, siginfo_t*, void*);

 protected:
  Scheduler() :
    cumulative_power_(0),
    num_power_samples_(0),
    max_power_(0),
    power_limit_(1000.0),
    available_power_(power_limit_)
  {
    if (global){
      fprintf(stderr, "only allowed one instance of a scheduler at a time\n");
      abort();
    }
    global = this;
    // initialize power measurement
    if(mic_open_device(&mic_device_, 0) != E_MIC_SUCCESS){
      std::cerr << "Error: Unable to open MIC" << std::endl;
      abort();
    }
  }
 private:
  virtual void runMaster(Task* root) = 0;

  static Scheduler* global;

  CPUList available_cores_;

  struct mic_device* mic_device_;

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

#endif

