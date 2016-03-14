
class Scheduler
{
 public:
  void runWorker();

  void terminateWorkers();
  
  virtual void init(int argc, char** argv);

  virtual void runMaster(Task* root) = 0;

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
 
 private:
  int nworkers_;
  int rank_;
  int nproc_;

};

class BasicScheduler : public Scheduler
{
 public:
  void runMaster(Task* root);
};

