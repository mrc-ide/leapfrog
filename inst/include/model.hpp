#pragma once

#include <consts.hpp>
#include <leapfrog.hpp>
#include <typedefs.hpp>
#include <utility>

template <typename real_type>
class Model {
 private:
  const Parameters pars;

  State state_cur;
  state state_next;

  template <typename real_type>
  void step() {
    run_demographic_projection<real_type>(pars, state_cur, state_next);
    std::swap(state_cur, state_next)
  }

 public:
  template <typename real_type>
  Model(Parameters& pars, State& initial_state;)
      : pars(pars), state_cur(initial_state) {
    // initialise population

    for (int g = 0; g < pars.num_genders; g++) {
      for (int a = 0; a < pars.age_groups_pop; a++) {
        state_cur.total_population(a, g) = pars.base_pop(a, g);
      }
    }
    state_next = state_cur;
  }

  template <typename real_type>
  void run_model(sim_years) {
    for (int time_step = 1; time_step < pars.sim_years; time_step++) {
      step();
      // TODO: report out at specified point e.g. only once every 10 steps?
      // always call report but reporter can manage how often it saves output
    }
  }
};