#pragma once

#include "intermediate_data.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void distribute_rate_over_sexes(const int t,
                                const Parameters<ModelVariant, real_type> &pars,
                                IntermediateData<ModelVariant, real_type> &intermediate) {
  const auto& p_in = pars.hiv.incidence;
  auto& i_ha = intermediate.hiv;

  real_type denominator = i_ha.hiv_neg_aggregate(MALE) +
                          p_in.relative_risk_sex(t) * i_ha.hiv_neg_aggregate(FEMALE);
  real_type total_neg = i_ha.hiv_neg_aggregate(MALE) + i_ha.hiv_neg_aggregate(FEMALE);
  i_ha.rate_sex(MALE) = p_in.total_rate(t) * total_neg / denominator;
  i_ha.rate_sex(FEMALE) = i_ha.rate_sex(MALE) * p_in.relative_risk_sex(t);
}

template<typename ModelVariant, typename real_type, bool OwnedData>
void run_calculate_incidence_rate(int t,
                                  const Parameters<ModelVariant, real_type> &pars,
                                  const State<ModelVariant, real_type, OwnedData> &state_curr,
                                  State<ModelVariant, real_type, OwnedData> &state_next,
                                  IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_op = pars.options;
  const auto& c_da = state_curr.dp;
  const auto& c_ha = state_curr.hiv;
  auto& i_ha = intermediate.hiv;

  const auto adult_incid_first_age_group = p_op.adult_incidence_first_age_group;
  const auto adult_incid_last_age_group = adult_incid_first_age_group + p_op.pAG_INCIDPOP;

  for (int g = 0; g < ss_d.NS; ++g) {
    for (int a = adult_incid_first_age_group; a < adult_incid_last_age_group; ++a) {
      i_ha.hiv_neg_aggregate(g) += c_da.p_total_pop(a, g) - c_ha.p_hiv_pop(a, g);
    }
  }

  distribute_rate_over_sexes<ModelVariant>(t, pars, intermediate);
}


template<typename ModelVariant, typename real_type, bool OwnedData>
void run_disease_progression_and_mortality(int hiv_step,
                                           int t,
                                           const Parameters<ModelVariant, real_type> &pars,
                                           const State<ModelVariant, real_type, OwnedData> &state_curr,
                                           State<ModelVariant, real_type, OwnedData> &state_next,
                                           IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_nh = pars.hiv.natural_history;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& i_ha = intermediate.hiv;

  for (int g = 0; g < ss_d.NS; ++g) {
    for (int ha = 0; ha < ss_h.hAG; ++ha) {
      for (int hm = 0; hm < ss_h.hDS; ++hm) {
        i_ha.cd4mx_scale = 1.0;
        if (p_nh.scale_cd4_mortality && t >= p_op.ts_art_start &&
            hm >= i_ha.everARTelig_idx && n_ha.h_hiv_adult(hm, ha, g) > 0.0) {
          i_ha.artpop_hahm = 0.0;
          for (int hu = 0; hu < ss_h.hTS; ++hu) {
            i_ha.artpop_hahm += n_ha.h_art_adult(hu, hm, ha, g);
          }
          i_ha.cd4mx_scale = n_ha.h_hiv_adult(hm, ha, g) /
                             (n_ha.h_hiv_adult(hm, ha, g) + i_ha.artpop_hahm);
        }

        i_ha.deaths = i_ha.cd4mx_scale * p_nh.cd4_mortality(hm, ha, g) * n_ha.h_hiv_adult(hm, ha, g);
        i_ha.p_hiv_deaths_age_sex(ha, g) += p_op.dt * i_ha.deaths;
        n_ha.h_hiv_deaths_no_art(hm, ha, g) += p_op.dt * i_ha.deaths;
        i_ha.grad(hm, ha, g) = -i_ha.deaths;
      }

      for (int hm = 1; hm < ss_h.hDS; ++hm) {
        const auto hiv_adults_progressing_cd4_stage = p_nh.cd4_progression(hm - 1, ha, g) * n_ha.h_hiv_adult(hm - 1, ha, g);
        i_ha.grad(hm - 1, ha, g) -= hiv_adults_progressing_cd4_stage;
        i_ha.grad(hm, ha, g) += hiv_adults_progressing_cd4_stage;
      }
    }
  }
}

template<typename ModelVariant, typename real_type, bool OwnedData>
void run_new_p_infections(int hiv_step,
                          int t,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type, OwnedData> &state_curr,
                          State<ModelVariant, real_type, OwnedData> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().dp;
  const auto& p_in = pars.hiv.incidence;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& n_dp = state_next.dp;
  auto& i_ha = intermediate.hiv;

  const auto adult_incid_first_age_group = p_op.adult_incidence_first_age_group;
  const auto adult_incid_last_age_group = adult_incid_first_age_group + p_op.pAG_INCIDPOP;

  for (int g = 0; g < ss.NS; ++g) {
    i_ha.hiv_negative_pop.setZero();
    i_ha.Xhivn_incagerr = 0.0;

    for (int a = adult_incid_first_age_group; a < ss.pAG; ++a) {
      i_ha.hiv_negative_pop(a) = n_dp.p_total_pop(a, g) - n_ha.p_hiv_pop(a, g);
    }

    for (int a = adult_incid_first_age_group; a < adult_incid_last_age_group; ++a) {
      i_ha.Xhivn_incagerr += p_in.relative_risk_age(a - adult_incid_first_age_group, g, t) *
                             i_ha.hiv_negative_pop(a);
    }

    for (int a = adult_incid_first_age_group; a < ss.pAG; ++a) {
      i_ha.p_infections_ts(a, g) = i_ha.hiv_negative_pop(a) *
                                   i_ha.rate_sex(g) *
                                   p_in.relative_risk_age(a - adult_incid_first_age_group, g, t) *
                                   i_ha.hiv_neg_aggregate(g) /
                                   i_ha.Xhivn_incagerr;
    }
  }
}

template<typename ModelVariant, typename real_type, bool OwnedData>
void run_new_hiv_p_infections(int hiv_step,
                              int t,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type, OwnedData> &state_curr,
                              State<ModelVariant, real_type, OwnedData> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_nh = pars.hiv.natural_history;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& i_ha = intermediate.hiv;

  for (int g = 0; g < ss_d.NS; g++) {
    int a = p_op.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss_h.hAG; ++ha) {
      i_ha.p_infections_ha = 0.0;
      for (int i = 0; i < ss_h.hAG_span[ha]; i++, a++) {
        i_ha.p_infections_a = i_ha.p_infections_ts(a, g);
        i_ha.p_infections_ha += i_ha.p_infections_a;
        const auto new_infections = p_op.dt * i_ha.p_infections_a;
        n_ha.p_infections(a, g) += new_infections;
        n_ha.p_hiv_pop(a, g) += new_infections;
      }

      // add p_infections to grad hivpop
      for (int hm = 0; hm < ss_h.hDS; ++hm) {
        i_ha.grad(hm, ha, g) += i_ha.p_infections_ha * p_nh.cd4_initial_distribution(hm, ha, g);
      }
    }
  }
}

template<typename ModelVariant, typename real_type, bool OwnedData>
void run_art_progression_and_mortality(int hiv_step,
                                       int t,
                                       const Parameters<ModelVariant, real_type> &pars,
                                       const State<ModelVariant, real_type, OwnedData> &state_curr,
                                       State<ModelVariant, real_type, OwnedData> &state_next,
                                       IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_tx = pars.hiv.art;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& i_ha = intermediate.hiv;

  for (int g = 0; g < ss_d.NS; ++g) {
    for (int ha = 0; ha < ss_h.hAG; ++ha) {
      for (int hm = i_ha.everARTelig_idx; hm < ss_h.hDS; ++hm) {
        for (int hu = 0; hu < ss_h.hTS; ++hu) {
          i_ha.deaths_art = p_tx.mortality(hu, hm, ha, g) *
                            p_tx.mortality_time_rate_ratio(hu, t) *
                            n_ha.h_art_adult(hu, hm, ha, g);
          const auto new_hiv_deaths_art = p_op.dt * i_ha.deaths_art;
          i_ha.p_hiv_deaths_age_sex(ha, g) += new_hiv_deaths_art;
          n_ha.h_hiv_deaths_art(hu, hm, ha, g) += new_hiv_deaths_art;
          i_ha.gradART(hu, hm, ha, g) = -i_ha.deaths_art;
        }

        for (int hu = 0; hu < (ss_h.hTS - 1); ++hu) {
          const auto art_adults_progressing_treatment_stage = n_ha.h_art_adult(hu, hm, ha, g) / p_tx.h_art_stage_dur(hu);
          i_ha.gradART(hu, hm, ha, g) -= art_adults_progressing_treatment_stage;
          i_ha.gradART(hu + 1, hm, ha, g) += art_adults_progressing_treatment_stage;
        }

        // ART dropout
        if (p_tx.dropout_rate(t) > 0) {
          for (int hu = 0; hu < ss_h.hTS; ++hu) {
            const auto art_adult_dropout = p_tx.dropout_rate(t) * n_ha.h_art_adult(hu, hm, ha, g);
	          if (p_tx.dropout_recover_cd4 && hu >= 2 && hm >= 1) {
	            // recover people on ART >1 year to one higher CD4 category
	            i_ha.grad(hm - 1, ha, g) += art_adult_dropout;
	          } else {
	            i_ha.grad(hm, ha, g) += art_adult_dropout;
	          }
            i_ha.gradART(hu, hm, ha, g) -= art_adult_dropout;
          }
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type, bool OwnedData>
void run_h_art_initiation(int hiv_step,
                          int t,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type, OwnedData> &state_curr,
                          State<ModelVariant, real_type, OwnedData> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_nh = pars.hiv.natural_history;
  const auto& p_tx = pars.hiv.art;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& i_ha = intermediate.hiv;

  for (int g = 0; g < ss_d.NS; ++g) {
    i_ha.Xart_15plus = 0.0;

    i_ha.artelig_hm.setZero();
    i_ha.Xartelig_15plus = 0.0;

    i_ha.expect_mort_artelig_hm.setZero();
    i_ha.expect_mort_artelig15plus = 0.0;

    for (int ha = p_op.hIDX_15PLUS; ha < ss_h.hAG; ++ha) {
      for (int hm = i_ha.everARTelig_idx; hm < ss_h.hDS; ++hm) {
        if (hm >= i_ha.anyelig_idx) {
          // TODO: Implement special population ART eligibility
          real_type prop_elig = 1.0;
	        real_type tmp_artelig = prop_elig * n_ha.h_hiv_adult(hm, ha, g);
          i_ha.artelig_hahm(hm, ha - p_op.hIDX_15PLUS) = tmp_artelig;
          i_ha.artelig_hm(hm) += tmp_artelig;
          i_ha.Xartelig_15plus += tmp_artelig;

	        real_type tmp_expect_mort = p_nh.cd4_mortality(hm, ha, g) * i_ha.artelig_hahm(hm, ha - p_op.hIDX_15PLUS);
	        i_ha.expect_mort_artelig_hm(hm) += tmp_expect_mort;
          i_ha.expect_mort_artelig15plus += tmp_expect_mort;
        }

        for (int hu = 0; hu < ss_h.hTS; ++hu) {
          i_ha.Xart_15plus += n_ha.h_art_adult(hu, hm, ha, g) +
                              p_op.dt * i_ha.gradART(hu, hm, ha, g);
	      }
      }
    }

    // calculate number on ART at end of ts, based on number or percent
    real_type art_interp_w = p_op.dt * (hiv_step + 1.0);
    if (p_op.proj_period_int == internal::PROJPERIOD_MIDYEAR && art_interp_w < 0.5) {
      if (!p_tx.adults_on_art_is_percent(g, t - 2) && !p_tx.adults_on_art_is_percent(g, t - 1)) {
        // case: both values are numbers
        i_ha.artnum_hts = (0.5 - art_interp_w) * p_tx.adults_on_art(g, t - 2) +
                          (art_interp_w + 0.5) * p_tx.adults_on_art(g, t - 1);
      } else if (p_tx.adults_on_art_is_percent(g, t - 2) && p_tx.adults_on_art_is_percent(g, t - 1)) {
        // case: both values are percentages
        i_ha.artcov_hts = (0.5 - art_interp_w) * p_tx.adults_on_art(g, t - 2) +
	                        (art_interp_w + 0.5) * p_tx.adults_on_art(g, t - 1);
        i_ha.artnum_hts = i_ha.artcov_hts * (i_ha.Xart_15plus + i_ha.Xartelig_15plus);
      } else if (!p_tx.adults_on_art_is_percent(g, t - 2) && p_tx.adults_on_art_is_percent(g, t - 1)) {
        // case: value is percentage only at time t - 1
        // transition from number to percentage
        i_ha.curr_coverage = i_ha.Xart_15plus / (i_ha.Xart_15plus + i_ha.Xartelig_15plus);
        i_ha.artcov_hts = i_ha.curr_coverage +
                          (p_tx.adults_on_art(g, t - 1) - i_ha.curr_coverage) *
                          p_op.dt / (0.5 - p_op.dt * hiv_step);
        // back to number
        i_ha.artnum_hts = i_ha.artcov_hts * (i_ha.Xart_15plus + i_ha.Xartelig_15plus);
      }
    } else {
      // If the projection period is calendar year (>= Spectrum v6.2),
      // this condition is always followed, and it interpolates between
      // end of last year and current year (+ 1.0).
      // If projection period was mid-year (<= Spectrum v6.19), the second
      // half of the projection year interpolates the first half of the
      // calendar year (e.g. hts 7/10 for 2019 interpolates December 2018
      // to December 2019)

      if (p_op.proj_period_int == internal::PROJPERIOD_MIDYEAR) {
	      art_interp_w -= 0.5;
      }

      if (!p_tx.adults_on_art_is_percent(g, t - 1) && !p_tx.adults_on_art_is_percent(g, t)) {
        // case: both values are numbers
        i_ha.artnum_hts = (1.0 - art_interp_w) * p_tx.adults_on_art(g, t - 1) +
	                        art_interp_w * p_tx.adults_on_art(g, t);
      } else if (p_tx.adults_on_art_is_percent(g, t - 1) && p_tx.adults_on_art_is_percent(g, t)) {
        // case: both values are percentages
        i_ha.artcov_hts = (1.0 - art_interp_w) * p_tx.adults_on_art(g, t - 1) +
	                        art_interp_w * p_tx.adults_on_art(g, t);
        // transition to number
        i_ha.artnum_hts = i_ha.artcov_hts * (i_ha.Xart_15plus + i_ha.Xartelig_15plus);
      } else if (!p_tx.adults_on_art_is_percent(g, t - 1) && p_tx.adults_on_art_is_percent(g, t)) {
        // case: value is percentage only at time t
        // transition from number to percentage
        i_ha.curr_coverage = i_ha.Xart_15plus / (i_ha.Xart_15plus + i_ha.Xartelig_15plus);
        i_ha.artcov_hts = i_ha.curr_coverage +
                          (p_tx.adults_on_art(g, t) - i_ha.curr_coverage) *
                          p_op.dt / (1.0 - art_interp_w + p_op.dt);
        // back to number
        i_ha.artnum_hts = i_ha.artcov_hts * (i_ha.Xart_15plus + i_ha.Xartelig_15plus);
      }
    }

    // Desired number to initiate on ART
    i_ha.artinit_hts = std::max(i_ha.artnum_hts - i_ha.Xart_15plus, 0.0);

    // Use mixture of eligibility and expected mortality for initiation distribution
    // Spectrum ART allocation is 2-step process
    // 1. Allocate by CD4 category (weighted by 'eligible' and 'expected mortality')
    // 2. Allocate by age groups (weighted only by eligibility)

    // Step 1: allocate ART by CD4 stage
    for (int hm = i_ha.anyelig_idx; hm < ss_h.hDS; ++hm) {
      auto eligibility_by_stage = (1.0 - p_tx.initiation_mortality_weight) *
	                                i_ha.artelig_hm(hm) /
                                  i_ha.Xartelig_15plus;
      auto expected_mortality_by_stage = p_tx.initiation_mortality_weight *
	                                       i_ha.expect_mort_artelig_hm(hm) /
                                         i_ha.expect_mort_artelig15plus;
      i_ha.artinit_hm(hm) = i_ha.artinit_hts * (eligibility_by_stage + expected_mortality_by_stage);
    }

    // Step 2: within CD4 category, allocate ART by age proportional to
    // eligibility
    for (int ha = p_op.hIDX_15PLUS; ha < ss_h.hAG; ++ha) {
      for (int hm = i_ha.anyelig_idx; hm < ss_h.hDS; ++hm) {
        if (i_ha.artelig_hm(hm) > 0.0) {
          i_ha.artinit_hahm = i_ha.artinit_hm(hm) *
	                            i_ha.artelig_hahm(hm, ha - p_op.hIDX_15PLUS) /
	                            i_ha.artelig_hm(hm);
          i_ha.artinit_hahm = std::min(i_ha.artinit_hahm, i_ha.artelig_hahm(hm, ha - p_op.hIDX_15PLUS));
          i_ha.artinit_hahm = std::min(i_ha.artinit_hahm,
                                       n_ha.h_hiv_adult(hm, ha, g) + p_op.dt * i_ha.grad(hm, ha, g));
          i_ha.grad(hm, ha, g) -= i_ha.artinit_hahm / p_op.dt;
          i_ha.gradART(ART0MOS, hm, ha, g) += i_ha.artinit_hahm / p_op.dt;
          n_ha.h_art_initiation(hm, ha, g) += i_ha.artinit_hahm;
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type, bool OwnedData>
void run_update_art_stratification(int hiv_step,
                                   int t,
                                   const Parameters<ModelVariant, real_type> &pars,
                                   const State<ModelVariant, real_type, OwnedData> &state_curr,
                                   State<ModelVariant, real_type, OwnedData> &state_next,
                                   IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& i_ha = intermediate.hiv;

  for (int g = 0; g < ss_d.NS; ++g) {
    for (int ha = 0; ha < ss_h.hAG; ++ha) {
      for (int hm = i_ha.everARTelig_idx; hm < ss_h.hDS; ++hm) {
        for (int hu = 0; hu < ss_h.hTS; ++hu) {
          n_ha.h_art_adult(hu, hm, ha, g) += p_op.dt * i_ha.gradART(hu, hm, ha, g);
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type, bool OwnedData>
void run_update_hiv_stratification(int hiv_step,
                                   int t,
                                   const Parameters<ModelVariant, real_type> &pars,
                                   const State<ModelVariant, real_type, OwnedData> &state_curr,
                                   State<ModelVariant, real_type, OwnedData> &state_next,
                                   IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& i_ha = intermediate.hiv;

  for (int g = 0; g < ss_d.NS; ++g) {
    for (int ha = 0; ha < ss_h.hAG; ++ha) {
      for (int hm = 0; hm < ss_h.hDS; ++hm) {
        n_ha.h_hiv_adult(hm, ha, g) += p_op.dt * i_ha.grad(hm, ha, g);
      }
    }
  }
}

template<typename ModelVariant, typename real_type, bool OwnedData>
void run_remove_p_hiv_deaths(int hiv_step,
                             int t,
                             const Parameters<ModelVariant, real_type> &pars,
                             const State<ModelVariant, real_type, OwnedData> &state_curr,
                             State<ModelVariant, real_type, OwnedData> &state_next,
                             IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& n_dp = state_next.dp;
  auto& i_ha = intermediate.hiv;

  for (int g = 0; g < ss_d.NS; ++g) {
    // sum HIV+ population size in each hivpop age group
    int a = p_op.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss_h.hAG; ++ha) {
      i_ha.hivpop_ha(ha) = 0.0;
      for (int i = 0; i < ss_h.hAG_span[ha]; ++i, ++a) {
        i_ha.hivpop_ha(ha) += n_ha.p_hiv_pop(a, g);
      }
    }

    // remove hivdeaths proportionally to age-distribution within each age group
    a = p_op.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss_h.hAG; ++ha) {
      if (i_ha.hivpop_ha(ha) > 0) {
        i_ha.hivqx_ha = i_ha.p_hiv_deaths_age_sex(ha, g) / i_ha.hivpop_ha(ha);
        for (int i = 0; i < ss_h.hAG_span[ha]; ++i, ++a) {
          i_ha.hivdeaths_a = n_ha.p_hiv_pop(a, g) * i_ha.hivqx_ha;
          n_ha.p_hiv_deaths(a, g) += i_ha.hivdeaths_a;
          n_dp.p_total_pop(a, g) -= i_ha.hivdeaths_a;
          n_ha.p_hiv_pop(a, g) -= i_ha.hivdeaths_a;
        }
      } else {
        a += ss_h.hAG_span[ha];
      }
    }
  }
}

}

template<typename ModelVariant, typename real_type, bool OwnedData>
void run_hiv_model_simulation(int t,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type, OwnedData> &state_curr,
                              State<ModelVariant, real_type, OwnedData> &state_next,
                              internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().hiv;
  const auto& p_tx = pars.hiv.art;
  const auto& p_op = pars.options;
  auto& i_ha = intermediate.hiv;

  i_ha.everARTelig_idx = p_tx.idx_hm_elig(t) < ss.hDS ? p_tx.idx_hm_elig(t) : ss.hDS;
  i_ha.anyelig_idx = p_tx.idx_hm_elig(t);

  internal::run_calculate_incidence_rate<ModelVariant>(
    t, pars, state_curr, state_next, intermediate
  );

  for (int hiv_step = 0; hiv_step < p_op.hts_per_year; ++hiv_step) {
    i_ha.grad.setZero();
    i_ha.gradART.setZero();
    i_ha.p_hiv_deaths_age_sex.setZero();

    internal::run_disease_progression_and_mortality<ModelVariant>(
      hiv_step, t, pars, state_curr, state_next, intermediate
    );

    internal::run_new_p_infections<ModelVariant>(
      hiv_step, t, pars, state_curr, state_next, intermediate
    );

    internal::run_new_hiv_p_infections<ModelVariant>(
      hiv_step, t, pars, state_curr, state_next, intermediate
    );

    if (t >= p_op.ts_art_start) {
      internal::run_art_progression_and_mortality<ModelVariant>(
        hiv_step, t, pars, state_curr, state_next, intermediate
      );

      internal::run_h_art_initiation<ModelVariant>(
        hiv_step, t, pars, state_curr, state_next, intermediate
      );

      internal::run_update_art_stratification<ModelVariant>(
        hiv_step, t, pars, state_curr, state_next, intermediate
      );
    }

    internal::run_update_hiv_stratification<ModelVariant>(
      hiv_step, t, pars, state_curr, state_next, intermediate
    );

    internal::run_remove_p_hiv_deaths<ModelVariant>(
      hiv_step, t, pars, state_curr, state_next, intermediate
    );
  }
}

}
