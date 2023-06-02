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

  intermediate.everARTelig_idx =
      pars.artcd4elig_idx(time_step) < pars.disease_stages ? pars.artcd4elig_idx(time_step) : pars.disease_stages;
  intermediate.anyelig_idx = pars.artcd4elig_idx(time_step);
  for (int hiv_step = 0; hiv_step < pars.hiv_steps_per_year; ++hiv_step) {
    run_disease_progression_and_mortality(time_step, pars, state_curr, state_next, intermediate, hiv_step);
    run_new_infections(time_step, pars, state_curr, state_next, intermediate, hiv_step);
    if (time_step > pars.time_art_start) {
      run_art_progression_and_mortality(time_step, pars, state_curr, state_next, intermediate, hiv_step);
      run_art_initiation(time_step, pars, state_curr, state_next, intermediate, hiv_step);
    }
  }
}

template<typename real_type>
void run_add_new_hiv_infections(int time_step,
                                const Parameters<real_type> &pars,
                                const State<real_type> &state_curr,
                                State<real_type> &state_next,
                                IntermediateData<real_type> &intermediate) {
  // TODO: Add EPP_DIRECTINCID switch

  // Calculating new infections once per year (like Spectrum)
  for (int g = 0; g < pars.num_genders; ++g) {
    for (int a = pars.adult_incidence_first_age_group; a < pars.age_groups_pop; a++) {
      intermediate.hiv_negative_pop(a, g) = state_curr.total_population(a, g) - state_curr.hiv_population(a, g);
    }
  }

  for (int g = 0; g < pars.num_genders; g++) {
    for (int a = pars.adult_incidence_first_age_group;
         a < pars.adult_incidence_first_age_group + pars.pAG_INCIDPOP; a++) {
      intermediate.hiv_neg_aggregate(g) += intermediate.hiv_negative_pop(a, g);
      intermediate.Xhivn_incagerr(g) +=
          pars.incidence_relative_risk_age(a - pars.adult_incidence_first_age_group, g, time_step) *
          intermediate.hiv_negative_pop(a, g);
    }
  }

  intermediate.incidence_rate_sex(MALE) =
      pars.incidence_rate(time_step) * (intermediate.hiv_neg_aggregate(MALE) + intermediate.hiv_neg_aggregate(FEMALE)) /
      (intermediate.hiv_neg_aggregate(MALE) +
       pars.incidence_relative_risk_sex(time_step) * intermediate.hiv_neg_aggregate(FEMALE));
  intermediate.incidence_rate_sex(FEMALE) =
      pars.incidence_rate(time_step) * pars.incidence_relative_risk_sex(time_step) *
      (intermediate.hiv_neg_aggregate(MALE) + intermediate.hiv_neg_aggregate(FEMALE)) /
      (intermediate.hiv_neg_aggregate(MALE) +
       pars.incidence_relative_risk_sex(time_step) * intermediate.hiv_neg_aggregate(FEMALE));

  for (int g = 0; g < pars.num_genders; g++) {
    for (int a = pars.hiv_adult_first_age_group; a < pars.age_groups_pop; a++) {
      intermediate.infections_ts(a, g) =
          intermediate.hiv_negative_pop(a, g) * intermediate.incidence_rate_sex(g) *
          pars.incidence_relative_risk_age(a - pars.adult_incidence_first_age_group, g, time_step) *
          intermediate.hiv_neg_aggregate(g) /
          intermediate.Xhivn_incagerr(g);
    }
  }
}

template<typename real_type>
void run_disease_progression_and_mortality(int time_step,
                                           const Parameters<real_type> &pars,
                                           const State<real_type> &state_curr,
                                           State<real_type> &state_next,
                                           IntermediateData<real_type> &intermediate,
                                           int hiv_step) {
  for (int g = 0; g < pars.num_genders; g++) {
    for (int ha = 0; ha < pars.age_groups_hiv; ha++) {
      for (int hm = 0; hm < pars.disease_stages; hm++) {
        // TODO: Mortality scaling not yet implemented
        if (pars.scale_cd4_mortality &
            (time_step >= pars.time_art_start) &
            (hm >= intermediate.everARTelig_idx) &
            (state_next.hiv_strat_adult(hm, ha, g) > 0.0)) {
          for (int hu = 0; hu < pars.treatment_stages; hu++) {
            intermediate.artpop_hahm += state_next.art_strat_adult(hu, hm, ha, g);
          }
          intermediate.cd4mx_scale = state_next.hiv_strat_adult(hm, ha, g) /
                                     (state_next.hiv_strat_adult(hm, ha, g) + intermediate.artpop_hahm);
        }

        intermediate.deaths =
            intermediate.cd4mx_scale * pars.cd4_mortality(hm, ha, g) * state_next.hiv_strat_adult(hm, ha, g);
        intermediate.hiv_deaths_age_sex(ha, g) += pars.dt * intermediate.deaths;
        state_next.aids_deaths_no_art(hm, ha, g) += pars.dt * intermediate.deaths;
        intermediate.grad(hm, ha, g) = -intermediate.deaths;
      }

      for (int hm = 1; hm < pars.disease_stages; hm++) {
        intermediate.grad(hm - 1, ha, g) -=
            pars.cd4_progression(hm - 1, ha, g) * state_next.hiv_strat_adult(hm - 1, ha, g);
        intermediate.grad(hm, ha, g) += pars.cd4_progression(hm - 1, ha, g) * state_next.hiv_strat_adult(hm - 1, ha, g);
      }
    }
  }
}

template<typename real_type>
void run_new_infections(int time_step,
                        const Parameters<real_type> &pars,
                        const State<real_type> &state_curr,
                        State<real_type> &state_next,
                        IntermediateData<real_type> &intermediate,
                        int hiv_step) {
  for (int g = 0; g < pars.num_genders; g++) {
    int a = pars.hiv_adult_first_age_group;
    for (int ha = 0; ha < pars.age_groups_hiv; ha++) {
      for (int i = 0; i < pars.hiv_age_groups_span(ha); i++, a++) {
        intermediate.infections_a = intermediate.infections_ts(a, g);
        intermediate.infections_ha += intermediate.infections_a;
        state_next.infections(a, g) += pars.dt * intermediate.infections_a;
        state_next.hiv_population(a, g) += pars.dt * intermediate.infections_a;
      }

      // add infections to grad hivpop
      for (int hm = 0; hm < pars.disease_stages; hm++) {
        intermediate.grad(hm, ha, g) += intermediate.infections_ha * pars.cd4_initdist(hm, ha, g);
      }
    }
  }
}

template<typename real_type>
void run_art_progression_and_mortality(int time_step,
                                       const Parameters<real_type> &pars,
                                       const State<real_type> &state_curr,
                                       State<real_type> &state_next,
                                       IntermediateData<real_type> &intermediate,
                                       int hiv_step) {
  for (int g = 0; g < pars.num_genders; g++) {
    for (int ha = 0; ha < pars.age_groups_hiv; ha++) {
      for (int hm = intermediate.everARTelig_idx; hm < pars.disease_stages; hm++) {

        for (int hu = 0; hu < pars.treatment_stages; hu++) {
          intermediate.deaths_art =
              pars.art_mortality(hu, hm, ha, g) * pars.artmx_timerr(hu, time_step) *
              state_next.art_strat_adult(hu, hm, ha, g);
          intermediate.hiv_deaths_age_sex(ha, g) += pars.dt * intermediate.deaths_art;
          state_next.aids_deaths_art(hu, hm, ha, g) += pars.dt * intermediate.deaths_art;
          intermediate.gradART(hu, hm, ha, g) = -intermediate.deaths_art;
        }

        for (int hu = 0; hu < (pars.treatment_stages - 1); hu++) {
          intermediate.gradART(hu, hm, ha, g) +=
              -state_next.art_strat_adult(hu, hm, ha, g) / pars.h_art_stage_dur(hu);
          intermediate.gradART(hu + 1, hm, ha, g) +=
              state_next.art_strat_adult(hu, hm, ha, g) / pars.h_art_stage_dur(hu);
        }

        // ART dropout
        if (pars.art_dropout(time_step) > 0) {
          for (int hu = 0; hu < pars.treatment_stages; hu++) {
            intermediate.grad(hm, ha, g) += pars.art_dropout(time_step) * state_next.art_strat_adult(hu, hm, ha, g);
            intermediate.gradART(hu, hm, ha, g) -=
                pars.art_dropout(time_step) * state_next.art_strat_adult(hu, hm, ha, g);
          }
        }
      }
    }
  }
}

template<typename real_type>
void run_art_initiation(int time_step,
                        const Parameters<real_type> &pars,
                        const State<real_type> &state_curr,
                        State<real_type> &state_next,
                        IntermediateData<real_type> &intermediate,
                        int hiv_step) {
  for (int g = 0; g < pars.num_genders; g++) {
    for (int ha = pars.hIDX_15PLUS; ha < pars.age_groups_hiv; ha++) {
      for (int hm = pars.everARTelig_idx; hm < pars.disease_stages; hm++) {
        if (hm >= intermediate.anyelig_idx) {
          // TODO: Implement special population ART eligibility
          Type prop_elig = 1.0;
          intermediate.Xartelig_15plus += intermediate.artelig_hahm(hm, ha - pars.hIDX_15PLUS) =
              prop_elig * state_curr.hiv_strat_adult(hm, ha, g);
          intermediate.expect_mort_artelig15plus +=
              pars.cd4_mortality(hm, ha, g) * intermediate.artelig_hahm(hm, ha - pars.hIDX_15PLUS);
        }
        for (int hu = 0; hu < pars.treatment_stages; hu++)
          intermediate.Xart_15plus +=
              state_curr.art_strat_adult(hu, hm, ha, g) + pars.dt * intermediate.gradART(hu, hm, ha, g);
      }
    } // loop over ha

    // calculate number on ART at end of ts, based on number or percent
    if (pars.dt * (hts + 1) < 0.5) {
      if ((!pars.art15plus_isperc(g, time_step - 2)) & (!pars.art15plus_isperc(g, time_step - 1))) {
        // Both values are numbers
        intermediate.artnum_hts =
            (0.5 - pars.dt * (hts + 1)) * pars.art15plus_num(g, t - 2) +
            (pars.dt * (hts + 1) + 0.5) * pars.art15plus_num(g, time_step - 1);
      } else if (pars.art15plus_isperc(g, time_step - 2) & pars.art15plus_isperc(g, time_step - 1)) {
        // Both values are percentages
        intermediate.artcov_hts =
            (0.5 - pars.dt * (hts + 1)) * pars.art15plus_num(g, time_step - 2) +
            (pars.dt * (hts + 1) + 0.5) * pars.art15plus_num(g, time_step - 1);
        intermediate.artnum_hts = intermediate.artcov_hts * (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
      } else if ((!pars.art15plus_isperc(g, time_step - 2)) & pars.art15plus_isperc(g, time_step - 1)) {
        // Transition from number to percentage
        intermediate.curr_coverage =
            intermediate.Xart_15plus / (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
        intermediate.artcov_hts = intermediate.curr_coverage +
                                  (pars.art15plus_num(g, time_step - 1) - intermediate.curr_coverage) * pars.dt /
                                  (0.5 - pars.dt * hts);
        intermediate.artnum_hts = intermediate.artcov_hts * (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
      }
    } else {
      if ((!pars.art15plus_isperc(g, time_step - 1)) & (!pars.art15plus_isperc(g, time_step))) {
        // Both values are numbers
        intermediate.artnum_hts =
            (1.5 - pars.dt * (hts + 1)) * pars.art15plus_num(g, time_step - 1) +
            (pars.dt * (hts + 1) - 0.5) * pars.art15plus_num(g, time_step);
      } else if (pars.art15plus_isperc(g, time_step - 1) & pars.art15plus_isperc(g, time_step)) {
        // Both values are percentages
        Type artcov_hts =
            (1.5 - pars.dt * (hts + 1)) * pars.art15plus_num(g, time_step - 1) +
            (pars.dt * (hts + 1) - 0.5) * pars.art15plus_num(g, time_step);
        intermediate.artnum_hts = artcov_hts * (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
      } else if ((!pars.art15plus_isperc(g, time_step - 1)) & pars.art15plus_isperc(g, time_step)) {
        // Transition from number to percentage
        intermediate.curr_coverage =
            intermediate.Xart_15plus / (intermediate.Xart_15plus + Xintermediate.artelig_15plus);
        intermediate.artcov_hts = intermediate.curr_coverage +
                                  (pars.art15plus_num(g, time_step) - intermediate.curr_coverage) * pars.dt /
                                  (1.5 - pars.dt * hts);
        intermediate.artnum_hts = intermediate.artcov_hts * (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
      }
    }

    // Desired number to initiate on ART
    Type artinit_hts =
        intermediate.artnum_hts > intermediate.Xart_15plus ? intermediate.artnum_hts - intermediate.Xart_15plus : 0.0;

    // Use mixture of eligibility and expected mortality for initiation distribution
    for (int ha = pars.hIDX_15PLUS; ha < pars.age_groups_hiv; ha++) {
      for (int hm = intermediate.anyelig_idx; hm < pars.disease_stages; hm++) {

        if (intermediate.Xartelig_15plus > 0.0) {
          Type artinit_hahm = artinit_hts * intermediate.artelig_hahm(hm, ha - pars.hIDX_15PLUS) *
                              ((1.0 - art_alloc_mxweight) / intermediate.Xartelig_15plus +
                               art_alloc_mxweight * pars.cd4_mortality(hm, ha, g) / expect_mort_artelig15plus);
          if (artinit_hahm > intermediate.artelig_hahm(hm, ha - pars.hIDX_15PLUS)) {
            artinit_hahm = aintermediate.rtelig_hahm(hm, ha - pars.hIDX_15PLUS);
          }
          if (artinit_hahm > hivstrat_adult(hm, ha, g, t) + dt * grad(hm, ha, g)) {
            artinit_hahm = hivstrat_adult(hm, ha, g, t) + dt * grad(hm, ha, g);
          }
          grad(hm, ha, g) -= artinit_hahm / dt;
          gradART(ART0MOS, hm, ha, g) += artinit_hahm / dt;
          artinit(hm, ha, g, t) += artinit_hahm;
        }
      }
    }
  }
}

}
}
