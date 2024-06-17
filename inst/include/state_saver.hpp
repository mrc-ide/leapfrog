#pragma once

#include "intermediate_data.hpp"
#include "model_variants.hpp"
#include "generated/state_saver_types.hpp"

namespace leapfrog {

template<typename ModelVariant, typename real_type>
struct OutputState {
  BaseModelOutputState<ModelVariant, real_type> base;
  ChildModelOutputState<ModelVariant, real_type> children;

  OutputState(int output_years) :
    base(output_years),
    children(output_years) {}
};

template<typename ModelVariant, typename real_type>
class StateSaver {
public:
  BaseModelStateSaver<ModelVariant, real_type> base;
  ChildModelStateSaver<ModelVariant, real_type> children;

  StateSaver(int time_steps,
             std::vector<int> save_steps) :
      save_steps(save_steps),
      full_state(save_steps.size()) {
    for (int step: save_steps) {
      if (step < 0) {
        std::stringstream ss;
        ss << "Output step must be at least 0, got '" << step << "'." << std::endl;
        throw std::runtime_error(ss.str());
      }
      if (step > time_steps) {
        std::stringstream ss;
        ss << "Output step can be at most number of time steps run which is '" << time_steps << "', got step '" << step
           << "'." << std::endl;
        throw std::runtime_error(ss.str());
      }
    }
  }


  void save_state(const State<ModelVariant, real_type, true> &state, int current_year) {
    for (size_t i = 0; i < save_steps.size(); ++i) {
      if (current_year == save_steps[i]) {
        base.save_state(full_state.base, i, state);
        children.save_state(full_state.children, i, state);
      }
    }
  }

  const OutputState<ModelVariant, real_type> &get_full_state() const {
    return full_state;
  }

private:
  std::vector<int> save_steps;
  OutputState<ModelVariant, real_type> full_state;
};

}
