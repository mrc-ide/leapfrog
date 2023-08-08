#pragma once

#include "types.hpp"
#include "model_variants.hpp"

namespace leapfrog {

template<typename ModelVariant, typename real_type>
struct BaseModelOutputState {
  Tensor3<real_type> p_total_pop;
  Tensor3<real_type> p_total_pop_natural_deaths;
  Tensor3<real_type> p_hiv_pop;
  Tensor3<real_type> p_hiv_pop_natural_deaths;
  Tensor4<real_type> h_hiv_adult;
  Tensor5<real_type> h_art_adult;
  Tensor1<real_type> births;
  Tensor4<real_type> h_hiv_deaths_no_art;
  Tensor3<real_type> p_infections;
  Tensor5<real_type> h_hiv_deaths_art;
  Tensor4<real_type> h_art_initiation;
  Tensor3<real_type> p_hiv_deaths;

  BaseModelOutputState(int no_output_years)
      : p_total_pop(StateSpace<ModelVariant>().base.pAG,
                    StateSpace<ModelVariant>().base.NS, no_output_years),
        p_total_pop_natural_deaths(StateSpace<ModelVariant>().base.pAG,
                                   StateSpace<ModelVariant>().base.NS, no_output_years),
        p_hiv_pop(StateSpace<ModelVariant>().base.pAG,
                  StateSpace<ModelVariant>().base.NS, no_output_years),
        p_hiv_pop_natural_deaths(StateSpace<ModelVariant>().base.pAG,
                                 StateSpace<ModelVariant>().base.NS, no_output_years),
        h_hiv_adult(StateSpace<ModelVariant>().base.hDS,
                    StateSpace<ModelVariant>().base.hAG,
                    StateSpace<ModelVariant>().base.NS,
                    no_output_years),
        h_art_adult(StateSpace<ModelVariant>().base.hTS,
                    StateSpace<ModelVariant>().base.hDS,
                    StateSpace<ModelVariant>().base.hAG,
                    StateSpace<ModelVariant>().base.NS,
                    no_output_years),
        births(no_output_years),
        h_hiv_deaths_no_art(StateSpace<ModelVariant>().base.hDS,
                            StateSpace<ModelVariant>().base.hAG,
                            StateSpace<ModelVariant>().base.NS, no_output_years),
        p_infections(StateSpace<ModelVariant>().base.pAG,
                     StateSpace<ModelVariant>().base.NS, no_output_years),
        h_hiv_deaths_art(StateSpace<ModelVariant>().base.hTS,
                         StateSpace<ModelVariant>().base.hDS,
                         StateSpace<ModelVariant>().base.hAG,
                         StateSpace<ModelVariant>().base.NS, no_output_years),
        h_art_initiation(StateSpace<ModelVariant>().base.hDS,
                         StateSpace<ModelVariant>().base.hAG,
                         StateSpace<ModelVariant>().base.NS, no_output_years),
        p_hiv_deaths(StateSpace<ModelVariant>().base.pAG,
                     StateSpace<ModelVariant>().base.NS, no_output_years) {
    p_total_pop.setZero();
    p_total_pop_natural_deaths.setZero();
    p_hiv_pop.setZero();
    p_hiv_pop_natural_deaths.setZero();
    h_hiv_adult.setZero();
    h_art_adult.setZero();
    births.setZero();
    h_hiv_deaths_no_art.setZero();
    p_infections.setZero();
    h_hiv_deaths_art.setZero();
    h_art_initiation.setZero();
    p_hiv_deaths.setZero();
  }
};

template<typename ModelVariant, typename real_type>
struct ChildModelOutputState {
  ChildModelOutputState(int no_output_years) {}
};

template<typename real_type>
struct ChildModelOutputState<ChildModel, real_type> {
  Tensor5<real_type> hc_hiv_pop;

  ChildModelOutputState(int no_output_years)
      : hc_hiv_pop(StateSpace<ChildModel>().base.hDS, StateSpace<ChildModel>().children.hTM,
                   StateSpace<ChildModel>().base.pAG,
                   StateSpace<ChildModel>().base.NS, no_output_years) {
    hc_hiv_pop.setZero();
  }
};

template<typename ModelVariant, typename real_type>
struct OutputState {
  BaseModelOutputState<ModelVariant, real_type> base;
  ChildModelOutputState<ModelVariant, real_type> children;

  OutputState(int no_output_years) :
      base(no_output_years),
      children(no_output_years) {}
};

template<typename ModelVariant, typename real_type>
class StateSaver {
public:
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


  void save_state(const State<ModelVariant, real_type> &state, int current_year) {
    for (size_t i = 0; i < save_steps.size(); ++i) {
      if (current_year == save_steps[i]) {
        full_state.base.p_total_pop.chip(i, full_state.base.p_total_pop.NumDimensions - 1) = state.base.p_total_pop;
        full_state.base.p_total_pop_natural_deaths.chip(i, full_state.base.p_total_pop_natural_deaths.NumDimensions -
                                                           1) = state.base.p_total_pop_natural_deaths;
        full_state.base.p_hiv_pop.chip(i, full_state.base.p_hiv_pop.NumDimensions - 1) = state.base.p_hiv_pop;
        full_state.base.p_hiv_pop_natural_deaths.chip(i, full_state.base.p_hiv_pop_natural_deaths.NumDimensions -
                                                         1) = state.base.p_hiv_pop_natural_deaths;
        full_state.base.h_hiv_adult.chip(i, full_state.base.h_hiv_adult.NumDimensions - 1) = state.base.h_hiv_adult;
        full_state.base.h_art_adult.chip(i, full_state.base.h_art_adult.NumDimensions - 1) = state.base.h_art_adult;
        full_state.base.births(i) = state.base.births;
        full_state.base.h_hiv_deaths_no_art.chip(i, full_state.base.h_hiv_deaths_no_art.NumDimensions -
                                                    1) = state.base.h_hiv_deaths_no_art;
        full_state.base.p_infections.chip(i, full_state.base.p_infections.NumDimensions - 1) = state.base.p_infections;
        full_state.base.h_hiv_deaths_art.chip(i, full_state.base.h_hiv_deaths_art.NumDimensions -
                                                 1) = state.base.h_hiv_deaths_art;
        full_state.base.h_art_initiation.chip(i, full_state.base.h_art_initiation.NumDimensions -
                                                 1) = state.base.h_art_initiation;
        full_state.base.p_hiv_deaths.chip(i, full_state.base.p_hiv_deaths.NumDimensions - 1) = state.base.p_hiv_deaths;
        if constexpr (ModelVariant::run_child_model) {
          full_state.children.hc_hiv_pop.chip(i, full_state.children.hc_hiv_pop.NumDimensions - 1) =
              state.children.hc_hiv_pop;
        }
        return;
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
