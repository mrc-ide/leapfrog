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
  for (int s = 0; s < ss.NS; ++s) {
    //less than 5 because there is a cd4 transition between ages 4 and 5
    for (int a = 1; a < ss.hc2_agestart; ++a) {
      for (int hd = 0; hd < ss.hc1DS; ++hd) {
        for (int cat = 0 ; cat < ss.hcTT; ++cat) {
          state_next.hc1_hiv_pop(hd, cat, a, s) += state_curr.hc1_hiv_pop(hd, cat, a-1, s) * demog.survival_probability(a, g, time_step);

        }
        for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.hc1_art_pop(dur, hd, a, s) += state_curr.hc1_art_pop(dur, hd, a-1, s) * demog.survival_probability(a, g, time_step);
        }

      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {
    for (int hd = 0; hd < ss.hc1DS; ++hd) {
      for (int hd_alt = 0; hd_alt < ss.hc2DS; ++hd_alt) {
        for (int cat = 0 ; cat < ss.hcTT; ++cat) {
          state_next.hc2_hiv_pop(hd_alt, cat, 0, s) +=  state_curr.hc1_hiv_pop(hd, cat, ss.hc1_ageend, s) * demog.survival_probability(ss.hc2_agestart, g, time_step) * cpars.hc_cd4_transition(hd_alt, hd);
        }
        for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.hc2_art_pop(dur, hd_alt, 0, s) += state_curr.hc1_art_pop(dur, hd, ss.hc1_ageend, s) * demog.survival_probability(ss.hc2_agestart, g, time_step) * cpars.hc_cd4_transition(hd_alt, hd);
        }
      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {
    for (int a = (ss.hc2_agestart + 1); a < pars.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss.hc2DS; ++hd) {
        for (int cat = 0 ; cat < ss.hcTT; ++cat) {
          state_next.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s) += state_curr.hc2_hiv_pop(hd, cat, a- ss.hc2_agestart-1, s) * demog.survival_probability(a, g, time_step);

        }
        for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.hc2_art_pop(dur, hd, a - ss.hc2_agestart, s) += state_curr.hc2_art_pop(dur, hd, a- ss.hc2_agestart-1, s) * demog.survival_probability(a, g, time_step);
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

  for (int s = 0; s < ss.NS; ++s) {
    // Run only first 5 age groups in total population 0, 1, 2, 3, 4
    for (int a = 0; a < ss.hc2_agestart; ++a) {
      if (cpars.hc_nosocomial(time_step) > 0) {
        state_next.p_infections(a, s) = cpars.hc_nosocomial(time_step) / (5.0 * ss.NS); //5.0 is used because we want to evenly distribute across the 5 age groups in 0-4
        state_next.p_hiv_pop(a, s) += state_next.p_infections(a, s);

        for (int hd = 0; hd < ss.hc1DS; ++hd) {
          // putting them all in perinatal hcTT to match spec nosocomial
          if (cpars.hc1_cd4_dist(hd) > 0) {
            state_next.hc1_hiv_pop(hd, 0, a, s) += state_next.p_infections(a, s) * cpars.hc1_cd4_dist(hd);
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


 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 0; hd < ss.hc1DS; ++hd) {
     for (int a = 0; a < ss.hc2_agestart; ++a) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_posthivmort(hd, cat, a, s) += state_next.hc1_hiv_pop(hd, cat, a, s) - (1.0 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc1_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a);
         //unsure whether this should be state_next or state_curr
        // state_next.aidsdeaths_noart_paed(hd, cat, a, s) += (1 - cpars.ctx_effect * cpars.ctx_val(t)) * state_curr.hc_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a)  ;
       }
     }
   }
 }


 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 0; hd < ss.hc2DS; ++hd) {
     for (int a = ss.hc2_agestart; a < pars.options.p_idx_fertility_first; ++a) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_posthivmort(hd, cat, a, s) += state_next.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s) - (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s) * cpars.hc2_cd4_mort(hd, cat, a - ss.hc2_agestart);
         //state_next.aidsdeaths_noart_paed(hd, cat, a, s) +=  (1 - cpars.ctx_effect * cpars.ctx_val(t)) * state_curr.hc_hiv_pop(hd, cat, a, s) * cpars.hc2_cd4_mort(hd, cat, a - 5); // output hiv deaths, aggregated across transmission category
       }
     }
   }
 }

 //progress through CD4 categories
 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 1; hd < ss.hc1DS; ++hd) {
     for (int a = 0; a < ss.hc2_agestart; ++a) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_grad(hd - 1, cat, a, s) -=  (intermediate.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1) + state_next.hc1_hiv_pop(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1)) / 2; //moving to next cd4 category
         intermediate.hc_grad(hd, cat, a, s) += (intermediate.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1) + state_next.hc1_hiv_pop(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1)) / 2; //moving into this cd4 category
       }
     }
   }
 }

 //progress through CD4 categories
 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 1; hd < ss.hc2DS; ++hd) {
     for (int a = ss.hc2_agestart; a < pars.options.p_idx_fertility_first; ++a) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_grad(hd - 1, cat, a, s) -= (intermediate.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc2_cd4_prog(hd - 1) + state_next.hc2_hiv_pop(hd - 1, cat, a - ss.hc2_agestart, s) * cpars.hc2_cd4_prog(hd - 1)) / 2; //moving to next cd4 category
         intermediate.hc_grad(hd, cat, a, s) += (intermediate.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc2_cd4_prog(hd - 1) + state_next.hc2_hiv_pop(hd - 1, cat, a - ss.hc2_agestart, s) * cpars.hc2_cd4_prog(hd - 1)) / 2; //moving into this cd4 category
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


  for (int s = 0; s < ss.NS; ++s) {
    for (int hd = 0; hd < ss.hc2DS; ++hd) {
      for (int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for (int cat = 0; cat < ss.hcTT; ++cat) {
          if(a < ss.hc2_agestart){
            intermediate.hc_grad(hd, cat, a, s) -= (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.hc1_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a)  ;
          }else{
            intermediate.hc_grad(hd, cat, a, s) -= (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s) * cpars.hc2_cd4_mort(hd, cat, a - ss.hc2_agestart)  ;

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
  for (int s = 0; s < ss.NS; ++s) {
    for (int hd = 0; hd < ss.hc1DS; ++hd) {
      for (int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for (int cat = 0; cat < ss.hcTT; ++cat) {
          if(a < ss.hc2_agestart){
            state_next.hc1_hiv_pop(hd, cat, a, s) +=intermediate.hc_grad(hd, cat, a, s);
          }else{
          state_next.hc2_hiv_pop(hd, cat, a- ss.hc2_agestart, s) +=intermediate.hc_grad(hd, cat, a, s);
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

