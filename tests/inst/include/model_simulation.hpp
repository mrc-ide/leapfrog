#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void distribute_rate_over_sexes(
    const int time_step,
    const Parameters<ModelVariant, real_type> &pars,
    IntermediateData<ModelVariant, real_type> &intermediate) {
  const auto incidence = pars.base.incidence;
  real_type denominator = intermediate.base.hiv_neg_aggregate(MALE) +
                          incidence.relative_risk_sex(time_step) *
                          intermediate.base.hiv_neg_aggregate(FEMALE);
  real_type total_neg = intermediate.base.hiv_neg_aggregate(MALE) +
                        intermediate.base.hiv_neg_aggregate(FEMALE);
  intermediate.base.rate_sex(MALE) =
      incidence.total_rate(time_step) * (total_neg) / denominator;
  intermediate.base.rate_sex(FEMALE) =
      incidence.total_rate(time_step) * incidence.relative_risk_sex(time_step) *
      (total_neg) / denominator;
}

template<typename ModelVariant, typename real_type>
void run_calculate_incidence_rate(int time_step,
                                  const Parameters<ModelVariant, real_type> &pars,
                                  const State<ModelVariant, real_type> &state_curr,
                                  State<ModelVariant, real_type> &state_next,
                                  IntermediateData<ModelVariant, real_type> &intermediate) {
  const auto adult_incidence_first_age_group = pars.base.options.adult_incidence_first_age_group;
  constexpr auto ss = StateSpace<ModelVariant>().base;

  for (int g = 0; g < ss.NS; ++g) {
    for (int a = adult_incidence_first_age_group;
         a < adult_incidence_first_age_group +
             pars.base.options.pAG_INCIDPOP; ++a) {
      intermediate.base.hiv_neg_aggregate(g) +=
          state_curr.base.p_total_pop(a, g) - state_curr.base.p_hiv_pop(a, g);
    }
  }

  distribute_rate_over_sexes<ModelVariant>(time_step, pars, intermediate);
}


template<typename ModelVariant, typename real_type>
void run_disease_progression_and_mortality(int hiv_step,
                                           int time_step,
                                           const Parameters<ModelVariant, real_type> &pars,
                                           const State<ModelVariant, real_type> &state_curr,
                                           State<ModelVariant, real_type> &state_next,
                                           IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  const auto natural_history = pars.base.natural_history;
  const auto dt = pars.base.options.dt;
  for (int g = 0; g < ss.NS; ++g) {
    for (int ha = 0; ha < ss.hAG; ++ha) {
      for (int hm = 0; hm < ss.hDS; ++hm) {
        intermediate.base.cd4mx_scale = 1.0;
        if (pars.base.natural_history.scale_cd4_mortality &&
            (time_step >= pars.base.options.ts_art_start) &&
            (hm >= intermediate.base.everARTelig_idx) &&
            (state_next.base.h_hiv_adult(hm, ha, g) > 0.0)) {
          intermediate.base.artpop_hahm = 0.0;
          for (int hu = 0; hu < ss.hTS; ++hu) {
            intermediate.base.artpop_hahm += state_next.base.h_art_adult(hu, hm,
                                                                         ha, g);
          }
          intermediate.base.cd4mx_scale =
              state_next.base.h_hiv_adult(hm, ha, g) /
              (state_next.base.h_hiv_adult(hm, ha, g) +
               intermediate.base.artpop_hahm);
        }

        intermediate.base.deaths =
            intermediate.base.cd4mx_scale *
            natural_history.cd4_mortality(hm, ha, g) *
            state_next.base.h_hiv_adult(hm, ha, g);
        intermediate.base.p_hiv_deaths_age_sex(ha, g) +=
            dt * intermediate.base.deaths;
        state_next.base.h_hiv_deaths_no_art(hm, ha, g) +=
            dt * intermediate.base.deaths;
        intermediate.base.grad(hm, ha, g) = -intermediate.base.deaths;
      }

      for (int hm = 1; hm < ss.hDS; ++hm) {
        intermediate.base.grad(hm - 1, ha, g) -=
            natural_history.cd4_progression(hm - 1, ha, g) *
            state_next.base.h_hiv_adult(hm - 1, ha, g);
        intermediate.base.grad(hm, ha, g) +=
            natural_history.cd4_progression(hm - 1, ha, g) *
            state_next.base.h_hiv_adult(hm - 1, ha, g);
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_new_p_infections(int hiv_step,
                          int time_step,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type> &state_curr,
                          State<ModelVariant, real_type> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  const auto incidence = pars.base.incidence;
  const auto adult_incidence_first_age_group = pars.base.options.adult_incidence_first_age_group;

  for (int g = 0; g < ss.NS; ++g) {
    intermediate.base.hiv_negative_pop.setZero();
    intermediate.base.Xhivn_incagerr = 0.0;

    for (int a = adult_incidence_first_age_group; a < ss.pAG; ++a) {
      intermediate.base.hiv_negative_pop(a) =
          state_next.base.p_total_pop(a, g) - state_next.base.p_hiv_pop(a, g);
    }

    for (int a = adult_incidence_first_age_group;
         a < adult_incidence_first_age_group +
             pars.base.options.pAG_INCIDPOP; ++a) {
      intermediate.base.Xhivn_incagerr +=
          incidence.relative_risk_age(a - adult_incidence_first_age_group, g,
                                      time_step) *
          intermediate.base.hiv_negative_pop(a);
    }

    for (int a = adult_incidence_first_age_group; a < ss.pAG; ++a) {
      intermediate.base.p_infections_ts(a, g) =
          intermediate.base.hiv_negative_pop(a) *
          intermediate.base.rate_sex(g) *
          incidence.relative_risk_age(a - adult_incidence_first_age_group, g,
                                      time_step) *
          intermediate.base.hiv_neg_aggregate(g) /
          intermediate.base.Xhivn_incagerr;

    }
  }
}

template<typename ModelVariant, typename real_type>
void run_new_hiv_p_infections(int hiv_step,
                              int time_step,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  const auto natural_history = pars.base.natural_history;

  for (int g = 0; g < ss.NS; g++) {
    int a = pars.base.options.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss.hAG; ++ha) {
      intermediate.base.p_infections_ha = 0.0;
      for (int i = 0; i < ss.hAG_span[ha]; i++, a++) {
        intermediate.base.p_infections_a = intermediate.base.p_infections_ts(a,
                                                                             g);
        intermediate.base.p_infections_ha += intermediate.base.p_infections_a;
        state_next.base.p_infections(a, g) +=
            pars.base.options.dt * intermediate.base.p_infections_a;
        state_next.base.p_hiv_pop(a, g) +=
            pars.base.options.dt * intermediate.base.p_infections_a;
      }


      // add p_infections to grad hivpop
      for (int hm = 0; hm < ss.hDS; ++hm) {
        intermediate.base.grad(hm, ha, g) +=
            intermediate.base.p_infections_ha *
            natural_history.cd4_initial_distribution(hm, ha, g);
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_art_progression_and_mortality(int hiv_step,
                                       int time_step,
                                       const Parameters<ModelVariant, real_type> &pars,
                                       const State<ModelVariant, real_type> &state_curr,
                                       State<ModelVariant, real_type> &state_next,
                                       IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  const auto art = pars.base.art;
  for (int g = 0; g < ss.NS; ++g) {
    for (int ha = 0; ha < ss.hAG; ++ha) {
      for (int hm = intermediate.base.everARTelig_idx; hm < ss.hDS; ++hm) {
        for (int hu = 0; hu < ss.hTS; ++hu) {
          intermediate.base.deaths_art =
              art.mortality(hu, hm, ha, g) *
              art.mortaility_time_rate_ratio(hu, time_step) *
              state_next.base.h_art_adult(hu, hm, ha, g);
          intermediate.base.p_hiv_deaths_age_sex(ha, g) +=
              pars.base.options.dt * intermediate.base.deaths_art;
          state_next.base.h_hiv_deaths_art(hu, hm, ha, g) +=
              pars.base.options.dt * intermediate.base.deaths_art;
          intermediate.base.gradART(hu, hm, ha,
                                    g) = -intermediate.base.deaths_art;
        }

        for (int hu = 0; hu < (ss.hTS - 1); ++hu) {
          intermediate.base.gradART(hu, hm, ha, g) +=
              -state_next.base.h_art_adult(hu, hm, ha, g) /
              art.h_art_stage_dur(hu);
          intermediate.base.gradART(hu + 1, hm, ha, g) +=
              state_next.base.h_art_adult(hu, hm, ha, g) /
              art.h_art_stage_dur(hu);
        }

        // ART dropout
        if (art.dropout(time_step) > 0) {
          for (int hu = 0; hu < ss.hTS; ++hu) {
            intermediate.base.grad(hm, ha, g) += art.dropout(time_step) *
                                                 state_next.base.h_art_adult(hu,
                                                                             hm,
                                                                             ha,
                                                                             g);
            intermediate.base.gradART(hu, hm, ha, g) -=
                art.dropout(time_step) *
                state_next.base.h_art_adult(hu, hm, ha, g);
          }
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_h_art_initiation(int hiv_step,
                          int time_step,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type> &state_curr,
                          State<ModelVariant, real_type> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  const auto natural_history = pars.base.natural_history;
  const auto art = pars.base.art;
  const auto dt = pars.base.options.dt;
  const auto hIDX_15PLUS = pars.base.options.hIDX_15PLUS;
  for (int g = 0; g < ss.NS; ++g) {
    intermediate.base.Xart_15plus = 0.0;
    intermediate.base.Xartelig_15plus = 0.0;
    intermediate.base.expect_mort_artelig15plus = 0.0;
    for (int ha = hIDX_15PLUS; ha < ss.hAG; ++ha) {
      for (int hm = intermediate.base.everARTelig_idx; hm < ss.hDS; ++hm) {
        if (hm >= intermediate.base.anyelig_idx) {
          // TODO: Implement special population ART eligibility
          real_type prop_elig = 1.0;
          intermediate.base.artelig_hahm(hm, ha - hIDX_15PLUS) =
              prop_elig * state_next.base.h_hiv_adult(hm, ha, g);
          intermediate.base.Xartelig_15plus += intermediate.base.artelig_hahm(
              hm, ha - hIDX_15PLUS);
          intermediate.base.expect_mort_artelig15plus +=
              natural_history.cd4_mortality(hm, ha, g) *
              intermediate.base.artelig_hahm(hm, ha - hIDX_15PLUS);
        }
        for (int hu = 0; hu < ss.hTS; ++hu)
          intermediate.base.Xart_15plus +=
              state_next.base.h_art_adult(hu, hm, ha, g) +
              dt * intermediate.base.gradART(hu, hm, ha, g);
      }
    }

    // calculate number on ART at end of ts, based on number or percent
    if (dt * (hiv_step + 1) < 0.5) {
      if (!art.adults_on_art_is_percent(g, time_step - 2) &&
          !art.adults_on_art_is_percent(g, time_step - 1)) {
        // Both values are numbers
        intermediate.base.artnum_hts =
            (0.5 - dt * (hiv_step + 1)) * art.adults_on_art(g, time_step - 2) +
            (dt * (hiv_step + 1) + 0.5) * art.adults_on_art(g, time_step - 1);
      } else if (art.adults_on_art_is_percent(g, time_step - 2) &&
                 art.adults_on_art_is_percent(g, time_step - 1)) {
        // Both values are percentages
        intermediate.base.artcov_hts =
            (0.5 - dt * (hiv_step + 1)) * art.adults_on_art(g, time_step - 2) +
            (dt * (hiv_step + 1) + 0.5) * art.adults_on_art(g, time_step - 1);
        intermediate.base.artnum_hts =
            intermediate.base.artcov_hts *
            (intermediate.base.Xart_15plus + intermediate.base.Xartelig_15plus);
      } else if (!art.adults_on_art_is_percent(g, time_step - 2) &&
                 art.adults_on_art_is_percent(g, time_step - 1)) {
        // Transition from number to percentage
        intermediate.base.curr_coverage =
            intermediate.base.Xart_15plus /
            (intermediate.base.Xart_15plus + intermediate.base.Xartelig_15plus);
        intermediate.base.artcov_hts = intermediate.base.curr_coverage +
                                       (art.adults_on_art(g, time_step - 1) -
                                        intermediate.base.curr_coverage) * dt /
                                       (0.5 - dt * hiv_step);
        intermediate.base.artnum_hts =
            intermediate.base.artcov_hts *
            (intermediate.base.Xart_15plus + intermediate.base.Xartelig_15plus);
      }
    } else {
      if (!art.adults_on_art_is_percent(g, time_step - 1) &&
          !art.adults_on_art_is_percent(g, time_step)) {
        // Both values are numbers
        intermediate.base.artnum_hts =
            (1.5 - dt * (hiv_step + 1)) * art.adults_on_art(g, time_step - 1) +
            (dt * (hiv_step + 1) - 0.5) * art.adults_on_art(g, time_step);
      } else if (art.adults_on_art_is_percent(g, time_step - 1) &&
                 art.adults_on_art_is_percent(g, time_step)) {
        // Both values are percentages
        intermediate.base.artcov_hts =
            (1.5 - dt * (hiv_step + 1)) * art.adults_on_art(g, time_step - 1) +
            (dt * (hiv_step + 1) - 0.5) * art.adults_on_art(g, time_step);
        intermediate.base.artnum_hts =
            intermediate.base.artcov_hts *
            (intermediate.base.Xart_15plus + intermediate.base.Xartelig_15plus);
      } else if (!art.adults_on_art_is_percent(g, time_step - 1) &&
                 art.adults_on_art_is_percent(g, time_step)) {
        // Transition from number to percentage
        intermediate.base.curr_coverage =
            intermediate.base.Xart_15plus /
            (intermediate.base.Xart_15plus + intermediate.base.Xartelig_15plus);
        intermediate.base.artcov_hts = intermediate.base.curr_coverage +
                                       (art.adults_on_art(g, time_step) -
                                        intermediate.base.curr_coverage) * dt /
                                       (1.5 - dt * hiv_step);
        intermediate.base.artnum_hts =
            intermediate.base.artcov_hts *
            (intermediate.base.Xart_15plus + intermediate.base.Xartelig_15plus);
      }
    }

    // Desired number to initiate on ART
    intermediate.base.artinit_hts =
        intermediate.base.artnum_hts > intermediate.base.Xart_15plus ?
        intermediate.base.artnum_hts -
        intermediate.base.Xart_15plus : 0.0;

    // Use mixture of eligibility and expected mortality for initiation distribution
    for (int ha = hIDX_15PLUS; ha < ss.hAG; ++ha) {
      for (int hm = intermediate.base.anyelig_idx; hm < ss.hDS; ++hm) {
        if (intermediate.base.Xartelig_15plus > 0.0) {
          intermediate.base.artinit_hahm =
              intermediate.base.artinit_hts *
              intermediate.base.artelig_hahm(hm, ha - hIDX_15PLUS) *
              ((1.0 - pars.base.art.initiation_mortality_weight) /
               intermediate.base.Xartelig_15plus +
               pars.base.art.initiation_mortality_weight *
               natural_history.cd4_mortality(hm, ha, g) /
               intermediate.base.expect_mort_artelig15plus);
          if (intermediate.base.artinit_hahm >
              intermediate.base.artelig_hahm(hm, ha - hIDX_15PLUS)) {
            intermediate.base.artinit_hahm = intermediate.base.artelig_hahm(hm,
                                                                            ha -
                                                                            hIDX_15PLUS);
          }
          if (intermediate.base.artinit_hahm >
              state_next.base.h_hiv_adult(hm, ha, g) +
              dt * intermediate.base.grad(hm, ha, g)) {
            intermediate.base.artinit_hahm =
                state_next.base.h_hiv_adult(hm, ha, g) +
                dt * intermediate.base.grad(hm, ha, g);
          }
          intermediate.base.grad(hm, ha, g) -=
              intermediate.base.artinit_hahm / dt;
          intermediate.base.gradART(ART0MOS, hm, ha, g) +=
              intermediate.base.artinit_hahm / dt;
          state_next.base.h_art_initiation(hm, ha,
                                           g) += intermediate.base.artinit_hahm;
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_update_art_stratification(int hiv_step,
                                   int time_step,
                                   const Parameters<ModelVariant, real_type> &pars,
                                   const State<ModelVariant, real_type> &state_curr,
                                   State<ModelVariant, real_type> &state_next,
                                   IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  for (int g = 0; g < ss.NS; ++g) {
    for (int ha = 0; ha < ss.hAG; ++ha) {
      for (int hm = intermediate.base.everARTelig_idx; hm < ss.hDS; ++hm) {
        for (int hu = 0; hu < ss.hTS; ++hu) {
          state_next.base.h_art_adult(hu, hm, ha, g) +=
              pars.base.options.dt * intermediate.base.gradART(hu, hm, ha, g);
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_update_hiv_stratification(int hiv_step,
                                   int time_step,
                                   const Parameters<ModelVariant, real_type> &pars,
                                   const State<ModelVariant, real_type> &state_curr,
                                   State<ModelVariant, real_type> &state_next,
                                   IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;

  for (int g = 0; g < ss.NS; ++g) {
    for (int ha = 0; ha < ss.hAG; ++ha) {
      for (int hm = 0; hm < ss.hDS; ++hm) {
        state_next.base.h_hiv_adult(hm, ha, g) +=
            pars.base.options.dt * intermediate.base.grad(hm, ha, g);
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_remove_p_hiv_deaths(int hiv_step,
                             int time_step,
                             const Parameters<ModelVariant, real_type> &pars,
                             const State<ModelVariant, real_type> &state_curr,
                             State<ModelVariant, real_type> &state_next,
                             IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  for (int g = 0; g < ss.NS; ++g) {
    // sum HIV+ population size in each hivpop age group
    int a = pars.base.options.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss.hAG; ++ha) {
      intermediate.base.hivpop_ha(ha) = 0.0;
      for (int i = 0; i < ss.hAG_span[ha]; ++i, ++a) {
        intermediate.base.hivpop_ha(ha) += state_next.base.p_hiv_pop(a, g);
      }
    }

    // remove hivdeaths proportionally to age-distribution within each age group
    a = pars.base.options.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss.hAG; ++ha) {
      if (intermediate.base.hivpop_ha(ha) > 0) {
        intermediate.base.hivqx_ha =
            intermediate.base.p_hiv_deaths_age_sex(ha, g) /
            intermediate.base.hivpop_ha(ha);
        for (int i = 0; i < ss.hAG_span[ha]; ++i, ++a) {
          intermediate.base.hivdeaths_a =
              state_next.base.p_hiv_pop(a, g) * intermediate.base.hivqx_ha;
          state_next.base.p_hiv_deaths(a, g) += intermediate.base.hivdeaths_a;
          state_next.base.p_total_pop(a, g) -= intermediate.base.hivdeaths_a;
          state_next.base.p_hiv_pop(a, g) -= intermediate.base.hivdeaths_a;
        }
      } else {
        a += ss.hAG_span[ha];
      }
    }
  }
}

}

template<typename ModelVariant, typename real_type>
void run_hiv_model_simulation(int time_step,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  const auto art = pars.base.art;

  intermediate.base.everARTelig_idx =
      art.idx_hm_elig(time_step) < ss.hDS ? art.idx_hm_elig(time_step) : ss.hDS;
  intermediate.base.anyelig_idx = art.idx_hm_elig(time_step);

  internal::run_calculate_incidence_rate<ModelVariant>(time_step, pars,
                                                       state_curr, state_next,
                                                       intermediate);


  for (int hiv_step = 0;
       hiv_step < pars.base.options.hts_per_year; ++hiv_step) {
    intermediate.base.grad.setZero();
    intermediate.base.gradART.setZero();
    intermediate.base.p_hiv_deaths_age_sex.setZero();
    internal::run_disease_progression_and_mortality<ModelVariant>(hiv_step,
                                                                  time_step,
                                                                  pars,
                                                                  state_curr,
                                                                  state_next,
                                                                  intermediate);

    internal::run_new_p_infections<ModelVariant>(hiv_step, time_step, pars,
                                                 state_curr, state_next,
                                                 intermediate);


    internal::run_new_hiv_p_infections<ModelVariant>(hiv_step, time_step, pars,
                                                     state_curr, state_next,
                                                     intermediate);

    if (time_step >= pars.base.options.ts_art_start) {
      internal::run_art_progression_and_mortality<ModelVariant>(hiv_step,
                                                                time_step, pars,
                                                                state_curr,
                                                                state_next,
                                                                intermediate);

      internal::run_h_art_initiation<ModelVariant>(hiv_step, time_step, pars,
                                                   state_curr, state_next,
                                                   intermediate);

      internal::run_update_art_stratification<ModelVariant>(hiv_step, time_step,
                                                            pars, state_curr,
                                                            state_next,
                                                            intermediate);

    }

    internal::run_update_hiv_stratification<ModelVariant>(hiv_step, time_step,
                                                          pars, state_curr,
                                                          state_next,
                                                          intermediate);

    internal::run_remove_p_hiv_deaths<ModelVariant>(hiv_step, time_step, pars,
                                                    state_curr, state_next,
                                                    intermediate);

  }
}

}
