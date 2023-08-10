#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void run_hiv_child_infections(int time_step,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");

  constexpr auto adult_ss = StateSpace<ModelVariant>().base;
  const auto children = pars.children.children;

  for (int g = 0; g < adult_ss.NS; ++g) {
    // Run only first 5 age groups in total population 0, 1, 2, 3, 4
    for (int af = 0; af < 5; ++af) {
      if (children.hc_nosocomial(time_step) > 0) {
        // Divide by 10 because we want to evenly distribute over 2 genders and 5 age groups
        state_next.base.p_infections(af, g) = children.hc_nosocomial(time_step) / (5 * adult_ss.NS);
        state_next.base.p_hiv_pop(af, g) += state_next.base.p_infections(af, g);

        for (int hm = 0; hm < adult_ss.hDS; ++hm) {
          // putting them all in perinatal hTM to match spec nosocomial
          if (children.hc1_cd4_dist(hm) > 0) {
            state_next.children.hc_hiv_pop(hm, 0, af, g) +=
                state_next.base.p_infections(af, g) * children.hc1_cd4_dist(hm);
          }
        }
      }
    }
  }
}
}

template<typename ModelVariant, typename real_type>
void run_child_model_simulation(int time_step,
                                const Parameters<ModelVariant, real_type> &pars,
                                const State<ModelVariant, real_type> &state_curr,
                                State<ModelVariant, real_type> &state_next,
                                internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  internal::run_hiv_child_infections(time_step, pars, state_curr, state_next, intermediate);
}

}
