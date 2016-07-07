#include <iostream>
#include <miclib.h>
#include <unistd.h>
#include <signal.h>
#include <cstdlib>

bool done = false;

void myhandler(int){
  done = true;
}

int main(){
  done = false;
  struct sigaction sa;
  sa.sa_handler = myhandler;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = 0;
  sigaction(SIGINT, &sa, NULL);

  struct mic_device* mic_device_;
  if(mic_open_device(&mic_device_, 0) != E_MIC_SUCCESS){
    std::cerr << "Error: Unable to open MIC" << std::endl;
    return -1;
  }
  while(!done){
    struct mic_power_util_info* pinfo;
    if(mic_get_power_utilization_info(mic_device_, &pinfo) != E_MIC_SUCCESS){
      std::cerr << "Error: Unable to read power utilization info" << std::endl;
      return -1;
    }
    uint32_t power;
    if(mic_get_inst_power_readings(pinfo, &power) != E_MIC_SUCCESS){
      std::cerr << "Error: Unable to read power utilization info" << std::endl;
      return -1;
    }
    std::cout << power / 1000000.0 << "\n";
    mic_free_power_utilization_info(pinfo);
    usleep(250000);
  }
  mic_close_device(mic_device_);
  return 0;
}
