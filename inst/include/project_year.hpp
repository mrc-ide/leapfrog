#pragma once

#include "generated/parameter_types.hpp"
#include "generated/state_types.hpp"
#include "general_demographic_projection.hpp"
#include "hiv_demographic_projection.hpp"
#include "model_simulation.hpp"
#include "child_model_simulation.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void project_year(int time_step,
                  const Parameters<ModelVariant, real_type> &pars,
                  const State<ModelVariant, real_type> &state_curr,
                  State<ModelVariant, real_type> &state_next,
                  IntermediateData<ModelVariant, real_type> &intermediate) {
  run_general_pop_demographic_projection<ModelVariant>(time_step, pars, state_curr, state_next,
                                                       intermediate);
  run_hiv_pop_demographic_projection<ModelVariant>(time_step, pars, state_curr, state_next,
                                                   intermediate);
  run_hiv_model_simulation<ModelVariant>(time_step, pars, state_curr, state_next, intermediate);
  if constexpr (ModelVariant::run_child_model) {
    run_child_model_simulation<ModelVariant>(time_step, pars, state_curr, state_next, intermediate);
  }
}

} // namespace internal
} // namespace leapfrog
