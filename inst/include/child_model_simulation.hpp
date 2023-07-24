#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {

template<HivAgeStratification S, typename real_type>
void run_paed_ageing(int time_step,
                     const Parameters<real_type> &pars,
                     const State<S, real_type> &state_curr,
                     State<S, real_type> &state_next,
                     IntermediateData<S, real_type> &intermediate) {
  const auto demog = pars.demography;
  const auto children = pars.children;
  
  constexpr auto ss = StateSpace<S>();
  for (int g = 0; g < ss.num_genders; ++g) {
    //less than 5 because there is a cd4 transition between ages 4 and 5
    for (int ha = 1; ha < 5; ++ha) {
      for (int hm = 0; hm < ss.hC1_disease_stages; ++hm) {
        for (int cat = 0 ; cat < ss.hTM; ++cat) {
          state_next.hc_hiv_pop(hm, cat, ha, g) += state_curr.hc_hiv_pop(hm, cat, ha-1, g) * demog.survival(ha, g, time_step);
          
        }
        for(int dur = 0; dur < ss.treatment_stages; ++dur){
          state_next.hc_art_pop(dur, hm, ha, g) += state_curr.hc_art_pop(dur, hm, ha-1, g) * demog.survival(ha, g, time_step);
        }
        
      }
    }
  }
  
  for(int g = 0; g < ss.num_genders; ++g){
    for(int hm = 0; hm < ss.hC1_disease_stages; ++hm){
      for(int hm_alt = 0; hm_alt < ss.hC2_disease_stages; ++hm_alt){
        for(int cat = 0 ; cat < ss.hTM; ++cat){
          state_next.hc_hiv_pop(hm_alt, cat, 5, g) +=  state_curr.hc_hiv_pop(hm, cat, 4, g) * demog.survival(5, g, time_step) * children.hc_cd4_transition(hm_alt, hm);
        }
        for(int dur = 0; dur < ss.treatment_stages; ++dur){
          state_next.hc_art_pop(dur, hm_alt, 5, g) += state_curr.hc_art_pop(dur, hm, 4, g) * demog.survival(5, g, time_step) * children.hc_cd4_transition(hm_alt, hm);
        }
      }
    }
  }
  
  for(int g = 0; g < ss.num_genders; ++g){
    //thought it was ss.fertility_first_age_group but maybe not? ASK ROB
    for(int ha = 6; ha < 15; ++ha) {
      for(int hm = 0; hm < ss.hC2_disease_stages; ++hm){
        for(int cat = 0 ; cat < ss.hTM; ++cat){
          state_next.hc_hiv_pop(hm, cat, ha, g) += state_curr.hc_hiv_pop(hm, cat, ha-1, g) * demog.survival(ha, g, time_step);
          
        }
        for(int dur = 0; dur < ss.treatment_stages; ++dur){
          state_next.hc_art_pop(dur, hm, ha, g) += state_curr.hc_art_pop(dur, hm, ha-1, g) * demog.survival(ha, g, time_step);
        }
        
      }
    }
  }
  
  
  
}

template<HivAgeStratification S, typename real_type>
void run_hiv_child_infections(int time_step,
                              const Parameters<real_type> &pars,
                              const State<S, real_type> &state_curr,
                              State<S, real_type> &state_next,
                              IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto children = pars.children;

  for (int g = 0; g < ss.num_genders; ++g) {
    // Run only first 5 age groups in total population 0, 1, 2, 3, 4
    for (int af = 0; af < 5; ++af) {
      if (children.hc_nosocomial(time_step) > 0) {
        // Divide by 10 because we want to evenly distribute over 2 genders and 5 age groups
        state_next.infections(af, g) = children.hc_nosocomial(time_step) / (5 * ss.num_genders);
        state_next.hiv_population(af, g) += state_next.infections(af, g);

        for (int hm = 0; hm < ss.disease_stages; ++hm) {
          // putting them all in perinatal hTM to match spec nosocomial
          if (children.hc1_cd4_dist(hm) > 0) {
            state_next.hc_hiv_pop(hm, 0, af, g) += state_next.infections(af, g) * children.hc1_cd4_dist(hm);
          }
        }
      }
    }
  }
}
}

template<HivAgeStratification S, typename real_type>
void run_child_model_simulation(int time_step,
                                const Parameters<real_type> &pars,
                                const State<S, real_type> &state_curr,
                                State<S, real_type> &state_next,
                                internal::IntermediateData<S, real_type> &intermediate) {
  internal::run_paed_ageing(time_step, pars, state_curr, state_next, intermediate);
  internal::run_hiv_child_infections(time_step, pars, state_curr, state_next, intermediate);
}

}
