#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {

template<HivAgeStratification S, typename real_type>
void run_child_ageing(int time_step,
                     const Parameters<real_type> &pars,
                     const State<S, real_type> &state_curr,
                     State<S, real_type> &state_next,
                     IntermediateData<S, real_type> &intermediate) {
  const auto demog = pars.demography;
  const auto cpars = pars.children;

  constexpr auto ss = StateSpace<S>();
  for (int g = 0; g < ss.NS; ++g) {
    //less than 5 because there is a cd4 transition between ages 4 and 5
    for (int af = 1; af < ss.hc2_agestart; ++af) {
      for (int hm = 0; hm < ss.hc1DS; ++hm) {
        for (int cat = 0 ; cat < ss.hcTT; ++cat) {
          state_next.hc1_hiv_pop(hm, cat, af, g) += state_curr.hc1_hiv_pop(hm, cat, af-1, g) * demog.survival_probability(af, g, time_step);

        }
        for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.hc1_art_pop(dur, hm, af, g) += state_curr.hc1_art_pop(dur, hm, af-1, g) * demog.survival_probability(af, g, time_step);
        }

      }
    }
  }

  for (int g = 0; g < ss.NS; ++g) {
    for (int hm = 0; hm < ss.hc1DS; ++hm) {
      for (int hm_alt = 0; hm_alt < ss.hc2DS; ++hm_alt) {
        for (int cat = 0 ; cat < ss.hcTT; ++cat) {
          state_next.hc2_hiv_pop(hm_alt, cat, 0, g) +=  state_curr.hc1_hiv_pop(hm, cat, ss.hc1_ageend, g) * demog.survival_probability(ss.hc2_agestart, g, time_step) * cpars.hc_cd4_transition(hm_alt, hm);
        }
        for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.hc2_art_pop(dur, hm_alt, 0, g) += state_curr.hc1_art_pop(dur, hm, ss.hc1_ageend, g) * demog.survival_probability(ss.hc2_agestart, g, time_step) * cpars.hc_cd4_transition(hm_alt, hm);
        }
      }
    }
  }

  for (int g = 0; g < ss.NS; ++g) {
    for (int af = (ss.hc2_agestart + 1); af < pars.options.p_idx_fertility_first; ++af) {
      for (int hm = 0; hm < ss.hc2DS; ++hm) {
        for (int cat = 0 ; cat < ss.hcTT; ++cat) {
          state_next.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g) += state_curr.hc2_hiv_pop(hm, cat, af- ss.hc2_agestart-1, g) * demog.survival_probability(af, g, time_step);

        }
        for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.hc2_art_pop(dur, hm, af - ss.hc2_agestart, g) += state_curr.hc2_art_pop(dur, hm, af- ss.hc2_agestart-1, g) * demog.survival_probability(af, g, time_step);
        }

      }
    }
  }



}

template<HivAgeStratification S, typename real_type>
void run_child_hiv_infections(int time_step,
                              const Parameters<real_type> &pars,
                              const State<S, real_type> &state_curr,
                              State<S, real_type> &state_next,
                              IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto cpars = pars.children;

  for (int g = 0; g < ss.NS; ++g) {
    // Run only first 5 age groups in total population 0, 1, 2, 3, 4
    for (int af = 0; af < ss.hc2_agestart; ++af) {
      if (cpars.hc_nosocomial(time_step) > 0) {
        state_next.p_infections(af, g) = cpars.hc_nosocomial(time_step) / (5.0 * ss.NS); //5.0 is used because we want to evenly distribute across the 5 age groups in 0-4
        state_next.p_hiv_pop(af, g) += state_next.p_infections(af, g);

        for (int hm = 0; hm < ss.hc1DS; ++hm) {
          // putting them all in perinatal hcTT to match spec nosocomial
          if (cpars.hc1_cd4_dist(hm) > 0) {
            state_next.hc1_hiv_pop(hm, 0, af, g) += state_next.p_infections(af, g) * cpars.hc1_cd4_dist(hm);
          }
        }
      }
    }
  }
}


template<HivAgeStratification S, typename real_type>
void run_child_natural_history(int time_step,
                              const Parameters<real_type> &pars,
                              const State<S, real_type> &state_curr,
                              State<S, real_type> &state_next,
                              IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto cpars = pars.children;


 for (int g = 0; g < ss.NS; ++g) {
   for (int hm = 0; hm < ss.hc1DS; ++hm) {
     for (int af = 0; af < ss.hc2_agestart; ++af) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_posthivmort(hm, cat, af, g) += state_next.hc1_hiv_pop(hm, cat, af, g) - (1.0 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc1_hiv_pop(hm, cat, af, g) * cpars.hc1_cd4_mort(hm, cat, af);
         //unsure whether this should be state_next or state_curr
        // state_next.aidsdeaths_noart_paed(hm, cat, af, g) += (1 - cpars.ctx_effect * cpars.ctx_val(t)) * state_curr.hc_hiv_pop(hm, cat, af, g) * cpars.hc1_cd4_mort(hm, cat, af)  ;
       }
     }
   }
 }


 for (int g = 0; g < ss.NS; ++g) {
   for (int hm = 0; hm < ss.hc2DS; ++hm) {
     for (int af = ss.hc2_agestart; af < pars.options.p_idx_fertility_first; ++af) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_posthivmort(hm, cat, af, g) += state_next.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g) - (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g) * cpars.hc2_cd4_mort(hm, cat, af - ss.hc2_agestart);
         //state_next.aidsdeaths_noart_paed(hm, cat, af, g) +=  (1 - cpars.ctx_effect * cpars.ctx_val(t)) * state_curr.hc_hiv_pop(hm, cat, af, g) * cpars.hc2_cd4_mort(hm, cat, af - 5); // output hiv deaths, aggregated across transmission category
       }
     }
   }
 }

 //progress through CD4 categories
 for (int g = 0; g < ss.NS; ++g) {
   for (int hm = 1; hm < ss.hc1DS; ++hm) {
     for (int af = 0; af < ss.hc2_agestart; ++af) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_grad(hm - 1, cat, af, g) -=  (intermediate.hc_posthivmort(hm - 1, cat, af, g) * cpars.hc1_cd4_prog(hm - 1) + state_next.hc1_hiv_pop(hm - 1, cat, af, g) * cpars.hc1_cd4_prog(hm - 1)) / 2; //moving to next cd4 category
         intermediate.hc_grad(hm, cat, af, g) += (intermediate.hc_posthivmort(hm - 1, cat, af, g) * cpars.hc1_cd4_prog(hm - 1) + state_next.hc1_hiv_pop(hm - 1, cat, af, g) * cpars.hc1_cd4_prog(hm - 1)) / 2; //moving into this cd4 category
       }
     }
   }
 }

 //progress through CD4 categories
 for (int g = 0; g < ss.NS; ++g) {
   for (int hm = 1; hm < ss.hc2DS; ++hm) {
     for (int af = ss.hc2_agestart; af < pars.options.p_idx_fertility_first; ++af) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_grad(hm - 1, cat, af, g) -= (intermediate.hc_posthivmort(hm - 1, cat, af, g) * cpars.hc2_cd4_prog(hm - 1) + state_next.hc2_hiv_pop(hm - 1, cat, af - ss.hc2_agestart, g) * cpars.hc2_cd4_prog(hm - 1)) / 2; //moving to next cd4 category
         intermediate.hc_grad(hm, cat, af, g) += (intermediate.hc_posthivmort(hm - 1, cat, af, g) * cpars.hc2_cd4_prog(hm - 1) + state_next.hc2_hiv_pop(hm - 1, cat, af - ss.hc2_agestart, g) * cpars.hc2_cd4_prog(hm - 1)) / 2; //moving into this cd4 category
       }
     }
   }
 }
}

template<HivAgeStratification S, typename real_type>
void run_child_hiv_mort(int time_step,
                               const Parameters<real_type> &pars,
                               const State<S, real_type> &state_curr,
                               State<S, real_type> &state_next,
                               IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto cpars = pars.children;


  for (int g = 0; g < ss.NS; ++g) {
    for (int hm = 0; hm < ss.hc2DS; ++hm) {
      for (int af = 0; af < pars.options.p_idx_fertility_first; ++af) {
        for (int cat = 0; cat < ss.hcTT; ++cat) {
          if(af < ss.hc2_agestart){
            intermediate.hc_grad(hm, cat, af, g) -= (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.hc1_hiv_pop(hm, cat, af, g) * cpars.hc1_cd4_mort(hm, cat, af)  ;
          }else{
            intermediate.hc_grad(hm, cat, af, g) -= (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g) * cpars.hc2_cd4_mort(hm, cat, af - ss.hc2_agestart)  ;

          }
        }
      }
    }
  }
}

template<HivAgeStratification S, typename real_type>
void add_child_grad(int time_step,
                         const Parameters<real_type> &pars,
                         const State<S, real_type> &state_curr,
                         State<S, real_type> &state_next,
                         IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto cpars = pars.children;

  //add on transitions
  for (int g = 0; g < ss.NS; ++g) {
    for (int hm = 0; hm < ss.hc1DS; ++hm) {
      for (int af = 0; af < pars.options.p_idx_fertility_first; ++af) {
        for (int cat = 0; cat < ss.hcTT; ++cat) {
          if(af < ss.hc2_agestart){
            state_next.hc1_hiv_pop(hm, cat, af, g) +=intermediate.hc_grad(hm, cat, af, g);
          }else{
          state_next.hc2_hiv_pop(hm, cat, af- ss.hc2_agestart, g) +=intermediate.hc_grad(hm, cat, af, g);
        }
      }
    }
  }

}

}
}// namespace internal

template<HivAgeStratification S, typename real_type>
void run_child_model_simulation(int time_step,
                                const Parameters<real_type> &pars,
                                const State<S, real_type> &state_curr,
                                State<S, real_type> &state_next,
                                internal::IntermediateData<S, real_type> &intermediate) {
  internal::run_child_ageing(time_step, pars, state_curr, state_next, intermediate);
  internal::run_child_hiv_infections(time_step, pars, state_curr, state_next, intermediate);
  internal::run_child_natural_history(time_step, pars, state_curr, state_next, intermediate);
  internal::run_child_hiv_mort(time_step, pars, state_curr, state_next, intermediate);
  internal::add_child_grad(time_step, pars, state_curr, state_next, intermediate);

}



} // namespace leapfrog

