#pragma once

class Model {
 private:
  int num_steps;
  virtual void step() = 0;

 public:
  virtual void run_model() = 0;
};

class BaseModel : public Model {
 private:
  void step() { run_demographic_projection() }

 public:
  BaseModel() {}

  void run_model() {
    for (int step = 0; step < num_steps; step++) {
      step()
      // report out at specified point i.e. only once every 10 steps?
    }
  }
}