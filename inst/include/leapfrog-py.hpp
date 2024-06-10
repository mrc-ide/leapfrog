#pragma once

#include "types.hpp"
#include "state_types.hpp"
#include "project_year.hpp"

namespace leapfrog {

/**
 * @brief Run a simulation model over a specified number of time steps.
 *
 * @tparam ModelVariant The variant of the model to be run, used for compile time switching.
 * @tparam real_type The data type used for real numbers in the simulation, usually a double.
 * @param time_step The time step to run, used in index parameters
 * @param pars The parameters required for running the simulation, read only.
 * @param state_curr The current state of the model, read only.
 * @param state_next The next state of the model.
 * @return None, updates state_next in place
 */
template<typename ModelVariant, typename real_type>
void project_single_year(int time_step,
                         const Parameters<ModelVariant, real_type> &pars,
                         const State<ModelVariant, real_type> &state_curr,
                         State<ModelVariant, real_type> &state_next) {

  internal::IntermediateData<ModelVariant, real_type> intermediate(pars.base.options.hAG_15plus);
  intermediate.reset();

  internal::project_year(time_step, pars, state_curr, state_next, intermediate);
}

} // namespace leapfrog
