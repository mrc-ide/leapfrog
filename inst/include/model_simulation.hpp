#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {

template<typename real_type>
void run_hiv_model_simulation(int time_step,
                              const Parameters<real_type> &pars,
                              const State<real_type> &state_curr,
                              State<real_type> &state_next,
                              IntermediateData<real_type> &intermediate) {
  run_add_new_hiv_infections(time_step, pars, state_curr, state_next, intermediate);
//  run_disease_progression(time_step, pars, state_curr, state_next, intermediate);
//  run_mortality(time_step, pars, state_curr, state_next, intermediate);
//  run_art_dropout(time_step, pars, state_curr, state_next, intermediate);
//  run_art_initiation(time_step, pars, state_curr, state_next, intermediate);
//  run_update_hiv_pop(time_step, pars, state_curr, state_next, intermediate);
//  run_hiv_deaths(time_step, pars, state_curr, state_next, intermediate);
}

template<typename real_type>
void run_add_new_hiv_infections(int time_step,
                                const Parameters<real_type> &pars,
                                const State<real_type> &state_curr,
                                State<real_type> &state_next,
                                IntermediateData<real_type> &intermediate) {
  // TODO: Add EPP_DIRECTINCID switch

  // Calculating new infections once per year (like Spectrum)
  Eigen::TensorFixedSize <real_type, Eigen::Sizes<pars.age_groups_pop, pars.num_genders>> hiv_negative_pop;
  for (int g = 0; g < pars.num_genders; ++g) {
    for (int a = pars.adult_incidence_first_age_group; a < pars.age_groups_pop; a++) {
      hiv_negative_pop(a, g) = state_curr.total_population(a, g) - state_curr.hiv_population(a, g);
    }
  }


  real_type hiv_neg_aggregate[pars.num_genders], Xhivn_incagerr[pars.num_genders];
  for (int g = 0; g < pars.num_genders; g++) {
    hiv_neg_aggregate[g] = 0.0;
    Xhivn_incagerr[g] = 0.0;
    for (int a = pars.adult_incidence_first_age_group;
         a < pars.adult_incidence_first_age_group + pars.pAG_INCIDPOP; a++) {
      hiv_neg_aggregate[g] += hiv_negative_pop(a, g);
      Xhivn_incagerr[g] +=
          pars.incidence_relative_risk_age(a - pars.adult_incidence_first_age_group, g, t) * hiv_negative_pop(a, g);
    }
  }

  real_type incidence_rate_sex[num_genders];
  incidence_rate_sex[MALE] = pars.incidence_rate(time_step) * (hiv_neg_aggregate[MALE] + hiv_neg_aggregate[FEMALE]) /
                             (hiv_neg_aggregate[MALE] +
                              pars.incidence_relative_risk_sex(time_step) * hiv_neg_aggregate[FEMALE]);
  incidence_rate_sex[FEMALE] =
      pars.incidence_rate(time_step) * pars.incidence_relative_risk_sex(time_step) *
      (hiv_neg_aggregate[MALE] + hiv_neg_aggregate[FEMALE]) /
      (hiv_neg_aggregate[MALE] + pars.incidence_relative_risk_sex(time_step) * hiv_neg_aggregate[FEMALE]);

  Eigen::TensorFixedSize <real_type, Eigen::Sizes<pars.age_groups_pop, pars.num_genders>> infections_ts;

  for (int g = 0; g < pars.num_genders; g++) {
    for (int a = pars.hiv_adult_first_age_group; a < pars.age_groups_pop; a++) {
      infections_ts(a, g) =
          hiv_negative_pop(a, g) * incidence_rate_sex[g] *
          pars.incidence_relative_risk_age(a - pars.adult_incidence_first_age_group, g, t) *
          hiv_neg_aggregate[g] /
          Xhivn_incagerr[g];
    }
  }
}

template<typename real_type>
void run_disease_progression_and_mortality(int time_step,
                                           const Parameters<real_type> &pars,
                                           const State<real_type> &state_curr,
                                           State<real_type> &state_next,
                                           IntermediateData<real_type> &intermediate) {

}


template<typename real_type>
void run_art_dropout(int time_step,
                     const Parameters<real_type> &pars,
                     const State<real_type> &state_curr,
                     State<real_type> &state_next,
                     IntermediateData<real_type> &intermediate) {

}

template<typename real_type>
void run_art_initiation(int time_step,
                        const Parameters<real_type> &pars,
                        const State<real_type> &state_curr,
                        State<real_type> &state_next,
                        IntermediateData<real_type> &intermediate) {

}

template<typename real_type>
void run_update_hiv_pop(int time_step,
                        const Parameters<real_type> &pars,
                        const State<real_type> &state_curr,
                        State<real_type> &state_next,
                        IntermediateData<real_type> &intermediate) {

}

template<typename real_type>
void run_hiv_deaths(int time_step,
                    const Parameters<real_type> &pars,
                    const State<real_type> &state_curr,
                    State<real_type> &state_next,
                    IntermediateData<real_type> &intermediate) {


}


}
}
