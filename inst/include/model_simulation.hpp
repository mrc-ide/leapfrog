#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {

template<typename real_type, HivAgeStratification S>
void distribute_incidence_rate_over_sexes(
    const int time_step,
    const Parameters<real_type> &pars,
    IntermediateData<real_type, S> &intermediate) {
  const auto incidence = pars.incidence;
  real_type denominator = intermediate.hiv_neg_aggregate(MALE) +
                          incidence.relative_risk_sex(time_step) * intermediate.hiv_neg_aggregate(FEMALE);
  real_type total_neg = intermediate.hiv_neg_aggregate(MALE) + intermediate.hiv_neg_aggregate(FEMALE);
  intermediate.incidence_rate_sex(MALE) = incidence.rate(time_step) * (total_neg) / denominator;
  intermediate.incidence_rate_sex(FEMALE) =
      incidence.rate(time_step) * incidence.relative_risk_sex(time_step) * (total_neg) / denominator;
}

template<typename real_type, HivAgeStratification S>
void run_add_new_hiv_infections(int time_step,
                                const Parameters<real_type> &pars,
                                const State<real_type, S> &state_curr,
                                State<real_type, S> &state_next,
                                IntermediateData<real_type, S> &intermediate) {
  // TODO: Add different HIV incidence rates see https://github.com/mrc-ide/leapfrog/issues/8

  const auto incidence = pars.incidence;
  const auto adult_incidence_first_age_group = pars.options.adult_incidence_first_age_group;
  constexpr auto ss = StateSpace<S>();

  // Calculating new infections once per year (like Spectrum)
  for (int g = 0; g < ss.num_genders; ++g) {
    for (int a = adult_incidence_first_age_group; a < ss.age_groups_pop; a++) {
      intermediate.hiv_negative_pop(a, g) = state_curr.total_population(a, g) - state_curr.hiv_population(a, g);
    }
  }

  for (int g = 0; g < ss.num_genders; g++) {
    for (int a = adult_incidence_first_age_group;
         a < adult_incidence_first_age_group + pars.options.pAG_INCIDPOP; a++) {
      intermediate.hiv_neg_aggregate(g) += intermediate.hiv_negative_pop(a, g);
      intermediate.Xhivn_incagerr(g) +=
          incidence.relative_risk_age(a - adult_incidence_first_age_group, g, time_step) *
          intermediate.hiv_negative_pop(a, g);
    }
  }

  distribute_incidence_rate_over_sexes<real_type, S>(time_step, pars, intermediate);

  for (int g = 0; g < ss.num_genders; g++) {
    for (int a = pars.options.hiv_adult_first_age_group; a < ss.age_groups_pop; a++) {
      intermediate.infections_ts(a, g) =
          intermediate.hiv_negative_pop(a, g) * intermediate.incidence_rate_sex(g) *
          incidence.relative_risk_age(a - adult_incidence_first_age_group, g, time_step) *
          intermediate.hiv_neg_aggregate(g) /
          intermediate.Xhivn_incagerr(g);
    }
  }
}


template<typename real_type, HivAgeStratification S>
void run_disease_progression_and_mortality(int hiv_step,
                                           int time_step,
                                           const Parameters<real_type> &pars,
                                           const State<real_type, S> &state_curr,
                                           State<real_type, S> &state_next,
                                           IntermediateData<real_type, S> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto natural_history = pars.natural_history;
  const auto dt = pars.options.dt;
  for (int g = 0; g < ss.num_genders; g++) {
    for (int ha = 0; ha < ss.age_groups_hiv; ++ha) {
      for (int hm = 0; hm < ss.disease_stages; ++hm) {
        intermediate.cd4mx_scale = 1.0;
        if (pars.options.scale_cd4_mortality &&
            (time_step >= pars.options.time_art_start) &&
            (hm >= intermediate.everARTelig_idx) &&
            (state_next.hiv_strat_adult(hm, ha, g) > 0.0)) {
          intermediate.artpop_hahm = 0.0;
          for (int hu = 0; hu < ss.treatment_stages; ++hu) {
            intermediate.artpop_hahm += state_next.art_strat_adult(hu, hm, ha, g);
          }
          intermediate.cd4mx_scale = state_next.hiv_strat_adult(hm, ha, g) /
                                     (state_next.hiv_strat_adult(hm, ha, g) + intermediate.artpop_hahm);
        }

        intermediate.deaths =
            intermediate.cd4mx_scale * natural_history.cd4_mortality(hm, ha, g) * state_next.hiv_strat_adult(hm, ha, g);
        intermediate.hiv_deaths_age_sex(ha, g) += dt * intermediate.deaths;
        state_next.aids_deaths_no_art(hm, ha, g) += dt * intermediate.deaths;
        intermediate.grad(hm, ha, g) = -intermediate.deaths;
      }

      for (int hm = 1; hm < ss.disease_stages; ++hm) {
        intermediate.grad(hm - 1, ha, g) -=
            natural_history.cd4_progression(hm - 1, ha, g) * state_next.hiv_strat_adult(hm - 1, ha, g);
        intermediate.grad(hm, ha, g) +=
            natural_history.cd4_progression(hm - 1, ha, g) * state_next.hiv_strat_adult(hm - 1, ha, g);
      }
    }
  }
}

template<typename real_type, HivAgeStratification S>
void run_new_infections(int hiv_step,
                        int time_step,
                        const Parameters<real_type> &pars,
                        const State<real_type, S> &state_curr,
                        State<real_type, S> &state_next,
                        IntermediateData<real_type, S> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto natural_history = pars.natural_history;
  for (int g = 0; g < ss.num_genders; g++) {
    int a = pars.options.hiv_adult_first_age_group;
    for (int ha = 0; ha < ss.age_groups_hiv; ++ha) {
      intermediate.infections_ha = 0.0;
      for (int i = 0; i < ss.hiv_age_groups_span[ha]; i++, a++) {
        intermediate.infections_a = intermediate.infections_ts(a, g);
        intermediate.infections_ha += intermediate.infections_a;
        state_next.infections(a, g) += pars.options.dt * intermediate.infections_a;
        state_next.hiv_population(a, g) += pars.options.dt * intermediate.infections_a;
      }

      // add infections to grad hivpop
      for (int hm = 0; hm < ss.disease_stages; ++hm) {
        intermediate.grad(hm, ha, g) += intermediate.infections_ha * natural_history.cd4_initdist(hm, ha, g);
      }
    }
  }
}

template<typename real_type, HivAgeStratification S>
void run_art_progression_and_mortality(int hiv_step,
                                       int time_step,
                                       const Parameters<real_type> &pars,
                                       const State<real_type, S> &state_curr,
                                       State<real_type, S> &state_next,
                                       IntermediateData<real_type, S> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto art = pars.art;
  for (int g = 0; g < ss.num_genders; g++) {
    for (int ha = 0; ha < ss.age_groups_hiv; ha++) {
      for (int hm = intermediate.everARTelig_idx; hm < ss.disease_stages; hm++) {
        for (int hu = 0; hu < ss.treatment_stages; hu++) {
          intermediate.deaths_art =
              art.mortality(hu, hm, ha, g) * art.artmx_timerr(hu, time_step) *
              state_next.art_strat_adult(hu, hm, ha, g);
          intermediate.hiv_deaths_age_sex(ha, g) += pars.options.dt * intermediate.deaths_art;
          state_next.aids_deaths_art(hu, hm, ha, g) += pars.options.dt * intermediate.deaths_art;
          intermediate.gradART(hu, hm, ha, g) = -intermediate.deaths_art;
        }

        for (int hu = 0; hu < (ss.treatment_stages - 1); hu++) {
          intermediate.gradART(hu, hm, ha, g) +=
              -state_next.art_strat_adult(hu, hm, ha, g) / art.h_art_stage_dur(hu);
          intermediate.gradART(hu + 1, hm, ha, g) +=
              state_next.art_strat_adult(hu, hm, ha, g) / art.h_art_stage_dur(hu);
        }

        // ART dropout
        if (art.dropout(time_step) > 0) {
          for (int hu = 0; hu < ss.treatment_stages; hu++) {
            intermediate.grad(hm, ha, g) += art.dropout(time_step) * state_next.art_strat_adult(hu, hm, ha, g);
            intermediate.gradART(hu, hm, ha, g) -=
                art.dropout(time_step) * state_next.art_strat_adult(hu, hm, ha, g);
          }
        }
      }
    }
  }
}

template<typename real_type, HivAgeStratification S>
void run_art_initiation(int hiv_step,
                        int time_step,
                        const Parameters<real_type> &pars,
                        const State<real_type, S> &state_curr,
                        State<real_type, S> &state_next,
                        IntermediateData<real_type, S> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto natural_history = pars.natural_history;
  const auto art = pars.art;
  const auto dt = pars.options.dt;
  const auto hIDX_15PLUS = pars.options.hIDX_15PLUS;
  for (int g = 0; g < ss.num_genders; ++g) {
    intermediate.Xart_15plus = 0.0;
    intermediate.Xartelig_15plus = 0.0;
    intermediate.expect_mort_artelig15plus = 0.0;
    for (int ha = hIDX_15PLUS; ha < ss.age_groups_hiv; ++ha) {
      for (int hm = intermediate.everARTelig_idx; hm < ss.disease_stages; ++hm) {
        if (hm >= intermediate.anyelig_idx) {
          // TODO: Implement special population ART eligibility
          real_type prop_elig = 1.0;
          intermediate.artelig_hahm(hm, ha - hIDX_15PLUS) =
              prop_elig * state_next.hiv_strat_adult(hm, ha, g);
          intermediate.Xartelig_15plus += intermediate.artelig_hahm(hm, ha - hIDX_15PLUS);
          intermediate.expect_mort_artelig15plus +=
              natural_history.cd4_mortality(hm, ha, g) * intermediate.artelig_hahm(hm, ha - hIDX_15PLUS);
        }
        for (int hu = 0; hu < ss.treatment_stages; ++hu)
          intermediate.Xart_15plus +=
              state_next.art_strat_adult(hu, hm, ha, g) + dt * intermediate.gradART(hu, hm, ha, g);
      }
    }

    // calculate number on ART at end of ts, based on number or percent
    if (dt * (hiv_step + 1) < 0.5) {
      if (!art.art15plus_isperc(g, time_step - 2) && !art.art15plus_isperc(g, time_step - 1)) {
        // Both values are numbers
        intermediate.artnum_hts =
            (0.5 - dt * (hiv_step + 1)) * art.art15plus_num(g, time_step - 2) +
            (dt * (hiv_step + 1) + 0.5) * art.art15plus_num(g, time_step - 1);
      } else if (art.art15plus_isperc(g, time_step - 2) && art.art15plus_isperc(g, time_step - 1)) {
        // Both values are percentages
        intermediate.artcov_hts =
            (0.5 - dt * (hiv_step + 1)) * art.art15plus_num(g, time_step - 2) +
            (dt * (hiv_step + 1) + 0.5) * art.art15plus_num(g, time_step - 1);
        intermediate.artnum_hts = intermediate.artcov_hts * (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
      } else if (!art.art15plus_isperc(g, time_step - 2) && art.art15plus_isperc(g, time_step - 1)) {
        // Transition from number to percentage
        intermediate.curr_coverage =
            intermediate.Xart_15plus / (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
        intermediate.artcov_hts = intermediate.curr_coverage +
                                  (art.art15plus_num(g, time_step - 1) - intermediate.curr_coverage) * dt /
                                  (0.5 - dt * hiv_step);
        intermediate.artnum_hts = intermediate.artcov_hts * (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
      }
    } else {
      if (!art.art15plus_isperc(g, time_step - 1) && !art.art15plus_isperc(g, time_step)) {
        // Both values are numbers
        intermediate.artnum_hts =
            (1.5 - dt * (hiv_step + 1)) * art.art15plus_num(g, time_step - 1) +
            (dt * (hiv_step + 1) - 0.5) * art.art15plus_num(g, time_step);
      } else if (art.art15plus_isperc(g, time_step - 1) && art.art15plus_isperc(g, time_step)) {
        // Both values are percentages
        intermediate.artcov_hts =
            (1.5 - dt * (hiv_step + 1)) * art.art15plus_num(g, time_step - 1) +
            (dt * (hiv_step + 1) - 0.5) * art.art15plus_num(g, time_step);
        intermediate.artnum_hts = intermediate.artcov_hts * (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
      } else if (!art.art15plus_isperc(g, time_step - 1) && art.art15plus_isperc(g, time_step)) {
        // Transition from number to percentage
        intermediate.curr_coverage =
            intermediate.Xart_15plus / (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
        intermediate.artcov_hts = intermediate.curr_coverage +
                                  (art.art15plus_num(g, time_step) - intermediate.curr_coverage) * dt /
                                  (1.5 - dt * hiv_step);
        intermediate.artnum_hts = intermediate.artcov_hts * (intermediate.Xart_15plus + intermediate.Xartelig_15plus);
      }
    }

    // Desired number to initiate on ART
    intermediate.artinit_hts =
        intermediate.artnum_hts > intermediate.Xart_15plus ? intermediate.artnum_hts - intermediate.Xart_15plus : 0.0;

    // Use mixture of eligibility and expected mortality for initiation distribution
    for (int ha = hIDX_15PLUS; ha < ss.age_groups_hiv; ++ha) {
      for (int hm = intermediate.anyelig_idx; hm < ss.disease_stages; ++hm) {
        if (intermediate.Xartelig_15plus > 0.0) {
          intermediate.artinit_hahm = intermediate.artinit_hts * intermediate.artelig_hahm(hm, ha - hIDX_15PLUS) *
                                      ((1.0 - pars.options.art_alloc_mxweight) / intermediate.Xartelig_15plus +
                                       pars.options.art_alloc_mxweight * natural_history.cd4_mortality(hm, ha, g) /
                                       intermediate.expect_mort_artelig15plus);
          if (intermediate.artinit_hahm > intermediate.artelig_hahm(hm, ha - hIDX_15PLUS)) {
            intermediate.artinit_hahm = intermediate.artelig_hahm(hm, ha - hIDX_15PLUS);
          }
          if (intermediate.artinit_hahm >
              state_next.hiv_strat_adult(hm, ha, g) + dt * intermediate.grad(hm, ha, g)) {
            intermediate.artinit_hahm =
                state_next.hiv_strat_adult(hm, ha, g) + dt * intermediate.grad(hm, ha, g);
          }
          intermediate.grad(hm, ha, g) -= intermediate.artinit_hahm / dt;
          intermediate.gradART(ART0MOS, hm, ha, g) += intermediate.artinit_hahm / dt;
          state_next.art_initiation(hm, ha, g) += intermediate.artinit_hahm;
        }
      }
    }
  }
}

template<typename real_type, HivAgeStratification S>
void run_update_art_stratification(int hiv_step,
                                   int time_step,
                                   const Parameters<real_type> &pars,
                                   const State<real_type, S> &state_curr,
                                   State<real_type, S> &state_next,
                                   IntermediateData<real_type, S> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  for (int g = 0; g < ss.num_genders; ++g) {
    for (int ha = 0; ha < ss.age_groups_hiv; ++ha) {
      for (int hm = intermediate.everARTelig_idx; hm < ss.disease_stages; ++hm) {
        for (int hu = 0; hu < ss.treatment_stages; ++hu) {
          state_next.art_strat_adult(hu, hm, ha, g) += pars.options.dt * intermediate.gradART(hu, hm, ha, g);
        }
      }
    }
  }
}

template<typename real_type, HivAgeStratification S>
void run_update_hiv_stratification(int hiv_step,
                                   int time_step,
                                   const Parameters<real_type> &pars,
                                   const State<real_type, S> &state_curr,
                                   State<real_type, S> &state_next,
                                   IntermediateData<real_type, S> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  for (int g = 0; g < ss.num_genders; ++g) {
    for (int ha = 0; ha < ss.age_groups_hiv; ++ha) {
      for (int hm = 0; hm < ss.disease_stages; ++hm) {
        state_next.hiv_strat_adult(hm, ha, g) += pars.options.dt * intermediate.grad(hm, ha, g);
      }
    }
  }
}

template<typename real_type, HivAgeStratification S>
void run_remove_hiv_deaths(int hiv_step,
                           int time_step,
                           const Parameters<real_type> &pars,
                           const State<real_type, S> &state_curr,
                           State<real_type, S> &state_next,
                           IntermediateData<real_type, S> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  for (int g = 0; g < ss.num_genders; ++g) {
    // sum HIV+ population size in each hivpop age group
    int a = pars.options.hiv_adult_first_age_group;
    for (int ha = 0; ha < ss.age_groups_hiv; ++ha) {
      intermediate.hivpop_ha(ha) = 0.0;
      for (int i = 0; i < ss.hiv_age_groups_span[ha]; ++i, ++a) {
        intermediate.hivpop_ha(ha) += state_next.hiv_population(a, g);
      }
    }

    // remove hivdeaths proportionally to age-distribution within each age group
    a = pars.options.hiv_adult_first_age_group;
    for (int ha = 0; ha < ss.age_groups_hiv; ++ha) {
      if (intermediate.hivpop_ha(ha) > 0) {
        intermediate.hivqx_ha = intermediate.hiv_deaths_age_sex(ha, g) / intermediate.hivpop_ha(ha);
        for (int i = 0; i < ss.hiv_age_groups_span[ha]; ++i, ++a) {
          intermediate.hivdeaths_a = state_next.hiv_population(a, g) * intermediate.hivqx_ha;
          state_next.hiv_deaths(a, g) += intermediate.hivdeaths_a;
          state_next.total_population(a, g) -= intermediate.hivdeaths_a;
          state_next.hiv_population(a, g) -= intermediate.hivdeaths_a;
        }
      } else {
        a += ss.hiv_age_groups_span[ha];
      }
    }
  }
}

}

template<typename real_type, HivAgeStratification S>
void run_hiv_model_simulation(int time_step,
                              const Parameters<real_type> &pars,
                              const State<real_type, S> &state_curr,
                              State<real_type, S> &state_next,
                              internal::IntermediateData<real_type, S> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto art = pars.art;

  internal::run_add_new_hiv_infections<real_type, S>(time_step, pars, state_curr, state_next, intermediate);

  intermediate.everARTelig_idx =
      art.artcd4elig_idx(time_step) < ss.disease_stages ? art.artcd4elig_idx(time_step) : ss.disease_stages;
  intermediate.anyelig_idx = art.artcd4elig_idx(time_step);

  for (int hiv_step = 0; hiv_step < pars.options.hiv_steps_per_year; ++hiv_step) {
    intermediate.grad.setZero();
    intermediate.gradART.setZero();
    intermediate.hiv_deaths_age_sex.setZero();

    internal::run_disease_progression_and_mortality<real_type, S>(hiv_step, time_step, pars, state_curr, state_next,
                                                                  intermediate);
    internal::run_new_infections<real_type, S>(hiv_step, time_step, pars, state_curr, state_next, intermediate);
    if (time_step >= pars.options.time_art_start) {
      internal::run_art_progression_and_mortality<real_type, S>(hiv_step, time_step, pars, state_curr, state_next,
                                                                intermediate);
      internal::run_art_initiation<real_type, S>(hiv_step, time_step, pars, state_curr, state_next, intermediate);
      internal::run_update_art_stratification<real_type, S>(hiv_step, time_step, pars, state_curr, state_next,
                                                            intermediate);
    }
    internal::run_update_hiv_stratification<real_type, S>(hiv_step, time_step, pars, state_curr, state_next,
                                                          intermediate);
    internal::run_remove_hiv_deaths<real_type, S>(hiv_step, time_step, pars, state_curr, state_next, intermediate);
  }
}

}
