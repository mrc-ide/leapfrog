#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void run_child_ageing(int time_step,
                      const Parameters<ModelVariant, real_type> &pars,
                      const State<ModelVariant, real_type> &state_curr,
                      State<ModelVariant, real_type> &state_next,
                      IntermediateData<ModelVariant, real_type> &intermediate) {
  const auto demog = pars.base.demography;
  const auto cpars = pars.children.children;
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;

  for (int s = 0; s < ss.NS; ++s) {
    //less than 5 because there is a cd4 transition between ages 4 and 5
    for (int a = 1; a < hc_ss.hc2_agestart; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
         for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          state_next.children.hc1_hiv_pop(hd, cat, a, s) +=
              state_curr.children.hc1_hiv_pop(hd, cat, a - 1, s) * demog.survival_probability(a, s, time_step);
        }
          for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.children.hc1_art_pop(dur, hd, a, s) +=
              state_curr.children.hc1_art_pop(dur, hd, a - 1, s) * demog.survival_probability(a, s, time_step);
        }
      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {
    for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
      for (int hd_alt = 0; hd_alt < hc_ss.hc2DS; ++hd_alt) {
        for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          state_next.children.hc2_hiv_pop(hd_alt, cat, 0, s) +=
              state_curr.children.hc1_hiv_pop(hd, cat, hc_ss.hc1_ageend, s) *
              demog.survival_probability(hc_ss.hc2_agestart, s, time_step) * cpars.hc_cd4_transition(hd_alt, hd);
        }
          for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.children.hc2_art_pop(dur, hd_alt, 0, s) +=
              state_curr.children.hc1_art_pop(dur, hd, hc_ss.hc1_ageend, s) *
              demog.survival_probability(hc_ss.hc2_agestart, s, time_step) * cpars.hc_cd4_transition(hd_alt, hd);
        }
      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {
    for (int a = (hc_ss.hc2_agestart + 1); a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc2DS; ++hd) {
       for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) +=
              state_curr.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart - 1, s) *
              demog.survival_probability(a, s, time_step);
        }
          for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.children.hc2_art_pop(dur, hd, a - hc_ss.hc2_agestart, s) +=
              state_curr.children.hc2_art_pop(dur, hd, a - hc_ss.hc2_agestart - 1, s) *
              demog.survival_probability(a, s, time_step);
        }

      }
    }
  }


}

template<typename ModelVariant, typename real_type>
void run_child_hiv_infections(int time_step,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;

  const auto cpars = pars.children.children;

  for (int s = 0; s < ss.NS; ++s) {
    // Run only first 5 age groups in total population 0, 1, 2, 3, 4
    for (int a = 0; a < hc_ss.hc2_agestart; ++a) {
      if (cpars.hc_nosocomial(time_step) > 0) {
        //5.0 is used because we want to evenly distribute across the 5 age groups in 0-4
        state_next.base.p_infections(a, s) = cpars.hc_nosocomial(time_step) / (5.0 * ss.NS);
        state_next.base.p_hiv_pop(a, s) += state_next.base.p_infections(a, s);

        for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          // putting them all in perinatal hcTT to match spec nosocomial
          if (cpars.hc1_cd4_dist(hd) > 0) {
            state_next.children.hc1_hiv_pop(hd, 0, a, s) += state_next.base.p_infections(a, s) * cpars.hc1_cd4_dist(hd);

          }
        }
      }
    }
  }
}


template<typename ModelVariant, typename real_type>
void run_child_natural_history(int time_step,
                               const Parameters<ModelVariant, real_type> &pars,
                               const State<ModelVariant, real_type> &state_curr,
                               State<ModelVariant, real_type> &state_next,
                               IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;

  for (int s = 0; s < ss.NS; ++s) {
    for (int a = 0; a < hc_ss.hc2_agestart; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          intermediate.children.hc_posthivmort(hd, cat, a, s) =
              state_next.children.hc1_hiv_pop(hd, cat, a, s) * (1 - cpars.hc1_cd4_mort(hd, cat, a));
          //intermediate.children.hc_posthivmort(hd, cat, a, s) = state_next.children.hc1_hiv_pop(hd, cat, a, s) - (1.0 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.children.hc1_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a);
        }
      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {
    for (int a = hc_ss.hc2_agestart; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
        for (int hd = 0; hd < hc_ss.hc2DS; ++hd) {
          intermediate.children.hc_posthivmort(hd, cat, a, s) +=
              state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) -
              (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) *
              state_curr.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) *
              cpars.hc2_cd4_mort(hd, cat, a - hc_ss.hc2_agestart);
        }
      }
    }
  }

  //progress through CD4 categories
  for (int s = 0; s < ss.NS; ++s) {
    for (int a = 0; a < hc_ss.hc2_agestart; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
        for (int hd = 1; hd < hc_ss.hc1DS; ++hd) {
          intermediate.children.hc_grad(hd - 1, cat, a, s) -=
              (intermediate.children.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1) +
               state_next.children.hc1_hiv_pop(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1)) /
              2; //moving to next cd4 category
          intermediate.children.hc_grad(hd, cat, a, s) +=
              (intermediate.children.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1) +
               state_next.children.hc1_hiv_pop(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1)) /
              2; //moving into this cd4 category


        }
      }
    }
  }
  //progress through CD4 categories
  for (int s = 0; s < ss.NS; ++s) {
    for (int a = hc_ss.hc2_agestart; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
        for (int hd = 1; hd < hc_ss.hc2DS; ++hd) {
          intermediate.children.hc_grad(hd - 1, cat, a, s) -=
              (intermediate.children.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc2_cd4_prog(hd - 1) +
               state_next.children.hc2_hiv_pop(hd - 1, cat, a - hc_ss.hc2_agestart, s) * cpars.hc2_cd4_prog(hd - 1)) /
              2; //moving to next cd4 category
          intermediate.children.hc_grad(hd, cat, a, s) +=
              (intermediate.children.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc2_cd4_prog(hd - 1) +
               state_next.children.hc2_hiv_pop(hd - 1, cat, a - hc_ss.hc2_agestart, s) * cpars.hc2_cd4_prog(hd - 1)) /
              2; //moving into this cd4 category
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_child_hiv_mort(int time_step,
                        const Parameters<ModelVariant, real_type> &pars,
                        const State<ModelVariant, real_type> &state_curr,
                        State<ModelVariant, real_type> &state_next,
                        IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;

  for (int s = 0; s < ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
        for (int hd = 0; hd < hc_ss.hc2DS; ++hd) {
          if (a < hc_ss.hc2_agestart) {
            //  intermediate.children.hc_grad(hd, cat, a, s) -= (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.children.hc1_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a);
            intermediate.children.hc_grad(hd, cat, a, s) -=
                state_next.children.hc1_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a);
          } else {
            // intermediate.children.hc_grad(hd, cat, a, s) -= (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.children.hc2_hiv_pop(hd, cat,a - hc_ss.hc2_agestart, s) *  cpars.hc2_cd4_mort(hd, cat, a - hc_ss.hc2_agestart);
            intermediate.children.hc_grad(hd, cat, a, s) -=
                state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) *
                cpars.hc2_cd4_mort(hd, cat, a - hc_ss.hc2_agestart);
          }
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void add_child_grad(int time_step,
                    const Parameters<ModelVariant, real_type> &pars,
                    const State<ModelVariant, real_type> &state_curr,
                    State<ModelVariant, real_type> &state_next,
                    IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;

  //add on transitions
  for (int s = 0; s < ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
        if (a < hc_ss.hc2_agestart) {
          for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
            state_next.children.hc1_hiv_pop(hd, cat, a, s) += intermediate.children.hc_grad(hd, cat, a, s);
          }// end hc1DS
        } else {
          for (int hd = 0; hd < hc_ss.hc2DS; ++hd) {
            state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) +=
                intermediate.children.hc_grad(hd, cat, a, s);
          }//end hc2DS
        }//end else
      }//end a
    }// end s
  }
}

}// namespace internal

template<typename ModelVariant, typename real_type>
void run_child_model_simulation(int time_step,
                                const Parameters<ModelVariant, real_type> &pars,
                                const State<ModelVariant, real_type> &state_curr,
                                State<ModelVariant, real_type> &state_next,
                                internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  internal::run_child_ageing(time_step, pars, state_curr, state_next, intermediate);
  internal::run_child_hiv_infections(time_step, pars, state_curr, state_next, intermediate);
  internal::run_child_natural_history(time_step, pars, state_curr, state_next, intermediate);
  internal::run_child_hiv_mort(time_step, pars, state_curr, state_next, intermediate);
  internal::add_child_grad(time_step, pars, state_curr, state_next, intermediate);
}

} // namespace leapfrog
