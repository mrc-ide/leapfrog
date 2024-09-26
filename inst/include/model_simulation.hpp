#pragma once

#include "intermediate_data.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void distribute_rate_over_sexes(const int t,
                                const Parameters<ModelVariant, real_type> &pars,
                                IntermediateData<ModelVariant, real_type> &intermediate) {
  const auto& p_in = pars.base.incidence;
  auto& i_ba = intermediate.base;

  real_type denominator = i_ba.hiv_neg_aggregate(MALE) +
                          p_in.relative_risk_sex(t) * i_ba.hiv_neg_aggregate(FEMALE);
  real_type total_neg = i_ba.hiv_neg_aggregate(MALE) + i_ba.hiv_neg_aggregate(FEMALE);
  i_ba.rate_sex(MALE) = p_in.total_rate(t) * total_neg / denominator;
  i_ba.rate_sex(FEMALE) = i_ba.rate_sex(MALE) * p_in.relative_risk_sex(t);
}

template<typename ModelVariant, typename real_type>
void run_calculate_incidence_rate(int t,
                                  const Parameters<ModelVariant, real_type> &pars,
                                  const State<ModelVariant, real_type> &state_curr,
                                  State<ModelVariant, real_type> &state_next,
                                  IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_op = pars.base.options;
  const auto& c_ba = state_curr.base;
  auto& i_ba = intermediate.base;

  const auto adult_incid_first_age_group = p_op.adult_incidence_first_age_group;
  const auto adult_incid_last_age_group = adult_incid_first_age_group + p_op.pAG_INCIDPOP;

  for (int g = 0; g < ss_b.NS; ++g) {
    for (int a = adult_incid_first_age_group; a < adult_incid_last_age_group; ++a) {
      i_ba.hiv_neg_aggregate(g) += c_ba.p_total_pop(a, g) - c_ba.p_hiv_pop(a, g);
    }
  }

  distribute_rate_over_sexes<ModelVariant>(t, pars, intermediate);
}


template<typename ModelVariant, typename real_type>
void run_disease_progression_and_mortality(int hiv_step,
                                           int t,
                                           const Parameters<ModelVariant, real_type> &pars,
                                           const State<ModelVariant, real_type> &state_curr,
                                           State<ModelVariant, real_type> &state_next,
                                           IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_nh = pars.base.natural_history;
  const auto& p_op = pars.base.options;
  auto& n_ba = state_next.base;
  auto& i_ba = intermediate.base;

  for (int g = 0; g < ss_b.NS; ++g) {
    for (int ha = 0; ha < ss_b.hAG; ++ha) {
      for (int hm = 0; hm < ss_b.hDS; ++hm) {
        i_ba.cd4mx_scale = 1.0;
        if (p_nh.scale_cd4_mortality && t >= p_op.ts_art_start &&
            hm >= i_ba.everARTelig_idx && n_ba.h_hiv_adult(hm, ha, g) > 0.0) {
          i_ba.artpop_hahm = 0.0;
          for (int hu = 0; hu < ss_b.hTS; ++hu) {
            i_ba.artpop_hahm += n_ba.h_art_adult(hu, hm, ha, g);
          }
          i_ba.cd4mx_scale = n_ba.h_hiv_adult(hm, ha, g) /
                             (n_ba.h_hiv_adult(hm, ha, g) + i_ba.artpop_hahm);
        }

        i_ba.deaths = i_ba.cd4mx_scale * p_nh.cd4_mortality(hm, ha, g) * n_ba.h_hiv_adult(hm, ha, g);
        i_ba.p_hiv_deaths_age_sex(ha, g) += p_op.dt * i_ba.deaths;
        n_ba.h_hiv_deaths_no_art(hm, ha, g) += p_op.dt * i_ba.deaths;
        i_ba.grad(hm, ha, g) = -i_ba.deaths;
      }

      for (int hm = 1; hm < ss_b.hDS; ++hm) {
        const auto hiv_adults_progressing_cd4_stage = p_nh.cd4_progression(hm - 1, ha, g) * n_ba.h_hiv_adult(hm - 1, ha, g);
        i_ba.grad(hm - 1, ha, g) -= hiv_adults_progressing_cd4_stage;
        i_ba.grad(hm, ha, g) += hiv_adults_progressing_cd4_stage;
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_new_p_infections(int hiv_step,
                          int t,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type> &state_curr,
                          State<ModelVariant, real_type> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  const auto& p_in = pars.base.incidence;
  const auto& p_op = pars.base.options;
  auto& n_ba = state_next.base;
  auto& i_ba = intermediate.base;

  const auto adult_incid_first_age_group = p_op.adult_incidence_first_age_group;
  const auto adult_incid_last_age_group = adult_incid_first_age_group + p_op.pAG_INCIDPOP;

  for (int g = 0; g < ss.NS; ++g) {
    i_ba.hiv_negative_pop.setZero();
    i_ba.Xhivn_incagerr = 0.0;

    for (int a = adult_incid_first_age_group; a < ss.pAG; ++a) {
      i_ba.hiv_negative_pop(a) = n_ba.p_total_pop(a, g) - n_ba.p_hiv_pop(a, g);
    }

    for (int a = adult_incid_first_age_group; a < adult_incid_last_age_group; ++a) {
      i_ba.Xhivn_incagerr += p_in.relative_risk_age(a - adult_incid_first_age_group, g, t) *
                             i_ba.hiv_negative_pop(a);
    }

    for (int a = adult_incid_first_age_group; a < ss.pAG; ++a) {
      i_ba.p_infections_ts(a, g) = i_ba.hiv_negative_pop(a) *
                                   i_ba.rate_sex(g) *
                                   p_in.relative_risk_age(a - adult_incid_first_age_group, g, t) *
                                   i_ba.hiv_neg_aggregate(g) /
                                   i_ba.Xhivn_incagerr;
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_new_hiv_p_infections(int hiv_step,
                              int t,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_nh = pars.base.natural_history;
  const auto& p_op = pars.base.options;
  auto& n_ba = state_next.base;
  auto& i_ba = intermediate.base;

  for (int g = 0; g < ss_b.NS; g++) {
    int a = p_op.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss_b.hAG; ++ha) {
      i_ba.p_infections_ha = 0.0;
      for (int i = 0; i < ss_b.hAG_span[ha]; i++, a++) {
        i_ba.p_infections_a = i_ba.p_infections_ts(a, g);
        i_ba.p_infections_ha += i_ba.p_infections_a;
        const auto new_infections = p_op.dt * i_ba.p_infections_a;
        n_ba.p_infections(a, g) += new_infections;
        n_ba.p_hiv_pop(a, g) += new_infections;
      }

      // add p_infections to grad hivpop
      for (int hm = 0; hm < ss_b.hDS; ++hm) {
        i_ba.grad(hm, ha, g) += i_ba.p_infections_ha * p_nh.cd4_initial_distribution(hm, ha, g);
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_art_progression_and_mortality(int hiv_step,
                                       int t,
                                       const Parameters<ModelVariant, real_type> &pars,
                                       const State<ModelVariant, real_type> &state_curr,
                                       State<ModelVariant, real_type> &state_next,
                                       IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_tx = pars.base.art;
  const auto& p_op = pars.base.options;
  auto& n_ba = state_next.base;
  auto& i_ba = intermediate.base;

  for (int g = 0; g < ss_b.NS; ++g) {
    for (int ha = 0; ha < ss_b.hAG; ++ha) {
      for (int hm = i_ba.everARTelig_idx; hm < ss_b.hDS; ++hm) {
        for (int hu = 0; hu < ss_b.hTS; ++hu) {
          i_ba.deaths_art = p_tx.mortality(hu, hm, ha, g) *
                            p_tx.mortaility_time_rate_ratio(hu, t) *
                            n_ba.h_art_adult(hu, hm, ha, g);
          const auto new_hiv_deaths_art = p_op.dt * i_ba.deaths_art;
          i_ba.p_hiv_deaths_age_sex(ha, g) += new_hiv_deaths_art;
          n_ba.h_hiv_deaths_art(hu, hm, ha, g) += new_hiv_deaths_art;
          i_ba.gradART(hu, hm, ha, g) = -i_ba.deaths_art;
        }

        for (int hu = 0; hu < (ss_b.hTS - 1); ++hu) {
          const auto art_adults_progressing_treatment_stage = n_ba.h_art_adult(hu, hm, ha, g) / p_tx.h_art_stage_dur(hu);
          i_ba.gradART(hu, hm, ha, g) -= art_adults_progressing_treatment_stage;
          i_ba.gradART(hu + 1, hm, ha, g) += art_adults_progressing_treatment_stage;
        }

        // ART dropout
        if (p_tx.dropout_rate(t) > 0) {
          for (int hu = 0; hu < ss_b.hTS; ++hu) {
            const auto art_adult_dropout = p_tx.dropout_rate(t) * n_ba.h_art_adult(hu, hm, ha, g);
	          if (p_tx.dropout_recover_cd4 && hu >= 2 && hm >= 1) {
	            // recover people on ART >1 year to one higher CD4 category
	            i_ba.grad(hm - 1, ha, g) += art_adult_dropout;
	          } else {
	            i_ba.grad(hm, ha, g) += art_adult_dropout;
	          }
            i_ba.gradART(hu, hm, ha, g) -= art_adult_dropout;
          }
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_h_art_initiation(int hiv_step,
                          int t,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type> &state_curr,
                          State<ModelVariant, real_type> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_nh = pars.base.natural_history;
  const auto& p_tx = pars.base.art;
  const auto& p_op = pars.base.options;
  auto& n_ba = state_next.base;
  auto& i_ba = intermediate.base;

  for (int g = 0; g < ss_b.NS; ++g) {
    i_ba.Xart_15plus = 0.0;

    i_ba.artelig_hm.setZero();
    i_ba.Xartelig_15plus = 0.0;

    i_ba.expect_mort_artelig_hm.setZero();
    i_ba.expect_mort_artelig15plus = 0.0;

    for (int ha = p_op.hIDX_15PLUS; ha < ss_b.hAG; ++ha) {
      for (int hm = i_ba.everARTelig_idx; hm < ss_b.hDS; ++hm) {
        if (hm >= i_ba.anyelig_idx) {
          // TODO: Implement special population ART eligibility
          real_type prop_elig = 1.0;
	        real_type tmp_artelig = prop_elig * n_ba.h_hiv_adult(hm, ha, g);
          i_ba.artelig_hahm(hm, ha - p_op.hIDX_15PLUS) = tmp_artelig;
          i_ba.artelig_hm(hm) += tmp_artelig;
          i_ba.Xartelig_15plus += tmp_artelig;

	        real_type tmp_expect_mort = p_nh.cd4_mortality(hm, ha, g) * i_ba.artelig_hahm(hm, ha - p_op.hIDX_15PLUS);
	        i_ba.expect_mort_artelig_hm(hm) += tmp_expect_mort;
          i_ba.expect_mort_artelig15plus += tmp_expect_mort;
        }
	
        for (int hu = 0; hu < ss_b.hTS; ++hu) {
          i_ba.Xart_15plus += n_ba.h_art_adult(hu, hm, ha, g) +
                              p_op.dt * i_ba.gradART(hu, hm, ha, g);
	      }
      }
    }

    // calculate number on ART at end of ts, based on number or percent
    real_type art_interp_w = p_op.dt * (hiv_step + 1.0);
    if (p_op.proj_period_int == internal::PROJPERIOD_MIDYEAR && art_interp_w < 0.5) {
      if (!p_tx.adults_on_art_is_percent(g, t - 2) && !p_tx.adults_on_art_is_percent(g, t - 1)) {
        // case: both values are numbers
        i_ba.artnum_hts = (0.5 - art_interp_w) * p_tx.adults_on_art(g, t - 2) +
                          (art_interp_w + 0.5) * p_tx.adults_on_art(g, t - 1);
      } else if (p_tx.adults_on_art_is_percent(g, t - 2) && p_tx.adults_on_art_is_percent(g, t - 1)) {
        // case: both values are percentages
        i_ba.artcov_hts = (0.5 - art_interp_w) * p_tx.adults_on_art(g, t - 2) +
	                        (art_interp_w + 0.5) * p_tx.adults_on_art(g, t - 1);
        i_ba.artnum_hts = i_ba.artcov_hts * (i_ba.Xart_15plus + i_ba.Xartelig_15plus);
      } else if (!p_tx.adults_on_art_is_percent(g, t - 2) && p_tx.adults_on_art_is_percent(g, t - 1)) {
        // case: value is percentage only at time t - 1
        // transition from number to percentage
        i_ba.curr_coverage = i_ba.Xart_15plus / (i_ba.Xart_15plus + i_ba.Xartelig_15plus);
        i_ba.artcov_hts = i_ba.curr_coverage +
                          (p_tx.adults_on_art(g, t - 1) - i_ba.curr_coverage) *
                          p_op.dt / (0.5 - p_op.dt * hiv_step);
        // back to number
        i_ba.artnum_hts = i_ba.artcov_hts * (i_ba.Xart_15plus + i_ba.Xartelig_15plus);
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
        i_ba.artnum_hts = (1.0 - art_interp_w) * p_tx.adults_on_art(g, t - 1) +
	                        art_interp_w * p_tx.adults_on_art(g, t);
      } else if (p_tx.adults_on_art_is_percent(g, t - 1) && p_tx.adults_on_art_is_percent(g, t)) {
        // case: both values are percentages
        i_ba.artcov_hts = (1.0 - art_interp_w) * p_tx.adults_on_art(g, t - 1) +
	                        art_interp_w * p_tx.adults_on_art(g, t);
        // transition to number
        i_ba.artnum_hts = i_ba.artcov_hts * (i_ba.Xart_15plus + i_ba.Xartelig_15plus);	
      } else if (!p_tx.adults_on_art_is_percent(g, t - 1) && p_tx.adults_on_art_is_percent(g, t)) {
        // case: value is percentage only at time t
        // transition from number to percentage
        i_ba.curr_coverage = i_ba.Xart_15plus / (i_ba.Xart_15plus + i_ba.Xartelig_15plus);
        i_ba.artcov_hts = i_ba.curr_coverage +
                          (p_tx.adults_on_art(g, t) - i_ba.curr_coverage) *
                          p_op.dt / (1.0 - art_interp_w + p_op.dt);
        // back to number
        i_ba.artnum_hts = i_ba.artcov_hts * (i_ba.Xart_15plus + i_ba.Xartelig_15plus);
      }
    }

    // Desired number to initiate on ART
    i_ba.artinit_hts = std::max(i_ba.artnum_hts - i_ba.Xart_15plus, 0.0);

    // Use mixture of eligibility and expected mortality for initiation distribution
    // Spectrum ART allocation is 2-step process
    // 1. Allocate by CD4 category (weighted by 'eligible' and 'expected mortality')
    // 2. Allocate by age groups (weighted only by eligibility)

    // Step 1: allocate ART by CD4 stage
    for (int hm = i_ba.anyelig_idx; hm < ss_b.hDS; ++hm) {
      auto eligibility_by_stage = (1.0 - p_tx.initiation_mortality_weight) *
	                                i_ba.artelig_hm(hm) /
                                  i_ba.Xartelig_15plus;
      auto expected_mortality_by_stage = p_tx.initiation_mortality_weight *
	                                       i_ba.expect_mort_artelig_hm(hm) /
                                         i_ba.expect_mort_artelig15plus;
      i_ba.artinit_hm(hm) = i_ba.artinit_hts * (eligibility_by_stage + expected_mortality_by_stage);
    }

    // Step 2: within CD4 category, allocate ART by age proportional to
    // eligibility
    for (int ha = p_op.hIDX_15PLUS; ha < ss_b.hAG; ++ha) {
      for (int hm = i_ba.anyelig_idx; hm < ss_b.hDS; ++hm) {
        if (i_ba.artelig_hm(hm) > 0.0) {
          i_ba.artinit_hahm = i_ba.artinit_hm(hm) *
	                            i_ba.artelig_hahm(hm, ha - p_op.hIDX_15PLUS) /
	                            i_ba.artelig_hm(hm);
          i_ba.artinit_hahm = std::min(i_ba.artinit_hahm, i_ba.artelig_hahm(hm, ha - p_op.hIDX_15PLUS));
          i_ba.artinit_hahm = std::min(i_ba.artinit_hahm,
                                       n_ba.h_hiv_adult(hm, ha, g) + p_op.dt * i_ba.grad(hm, ha, g));
          i_ba.grad(hm, ha, g) -= i_ba.artinit_hahm / p_op.dt;
          i_ba.gradART(ART0MOS, hm, ha, g) += i_ba.artinit_hahm / p_op.dt;
          n_ba.h_art_initiation(hm, ha, g) += i_ba.artinit_hahm;
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_update_art_stratification(int hiv_step,
                                   int t,
                                   const Parameters<ModelVariant, real_type> &pars,
                                   const State<ModelVariant, real_type> &state_curr,
                                   State<ModelVariant, real_type> &state_next,
                                   IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_op = pars.base.options;
  auto& n_ba = state_next.base;
  auto& i_ba = intermediate.base;

  for (int g = 0; g < ss_b.NS; ++g) {
    for (int ha = 0; ha < ss_b.hAG; ++ha) {
      for (int hm = i_ba.everARTelig_idx; hm < ss_b.hDS; ++hm) {
        for (int hu = 0; hu < ss_b.hTS; ++hu) {
          n_ba.h_art_adult(hu, hm, ha, g) += p_op.dt * i_ba.gradART(hu, hm, ha, g);
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_update_hiv_stratification(int hiv_step,
                                   int t,
                                   const Parameters<ModelVariant, real_type> &pars,
                                   const State<ModelVariant, real_type> &state_curr,
                                   State<ModelVariant, real_type> &state_next,
                                   IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_op = pars.base.options;
  auto& n_ba = state_next.base;
  auto& i_ba = intermediate.base;

  for (int g = 0; g < ss_b.NS; ++g) {
    for (int ha = 0; ha < ss_b.hAG; ++ha) {
      for (int hm = 0; hm < ss_b.hDS; ++hm) {
        n_ba.h_hiv_adult(hm, ha, g) += p_op.dt * i_ba.grad(hm, ha, g);
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_remove_p_hiv_deaths(int hiv_step,
                             int t,
                             const Parameters<ModelVariant, real_type> &pars,
                             const State<ModelVariant, real_type> &state_curr,
                             State<ModelVariant, real_type> &state_next,
                             IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_op = pars.base.options;
  auto& n_ba = state_next.base;
  auto& i_ba = intermediate.base;

  for (int g = 0; g < ss_b.NS; ++g) {
    // sum HIV+ population size in each hivpop age group
    int a = p_op.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss_b.hAG; ++ha) {
      i_ba.hivpop_ha(ha) = 0.0;
      for (int i = 0; i < ss_b.hAG_span[ha]; ++i, ++a) {
        i_ba.hivpop_ha(ha) += n_ba.p_hiv_pop(a, g);
      }
    }

    // remove hivdeaths proportionally to age-distribution within each age group
    a = p_op.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss_b.hAG; ++ha) {
      if (i_ba.hivpop_ha(ha) > 0) {
        i_ba.hivqx_ha = i_ba.p_hiv_deaths_age_sex(ha, g) / i_ba.hivpop_ha(ha);
        for (int i = 0; i < ss_b.hAG_span[ha]; ++i, ++a) {
          i_ba.hivdeaths_a = n_ba.p_hiv_pop(a, g) * i_ba.hivqx_ha;
          n_ba.p_hiv_deaths(a, g) += i_ba.hivdeaths_a;
          n_ba.p_total_pop(a, g) -= i_ba.hivdeaths_a;
          n_ba.p_hiv_pop(a, g) -= i_ba.hivdeaths_a;
        }
      } else {
        a += ss_b.hAG_span[ha];
      }
    }
  }
}

}

template<typename ModelVariant, typename real_type>
void run_hiv_model_simulation(int t,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  const auto& p_tx = pars.base.art;
  const auto& p_op = pars.base.options;
  auto& i_ba = intermediate.base;

  i_ba.everARTelig_idx = p_tx.idx_hm_elig(t) < ss.hDS ? p_tx.idx_hm_elig(t) : ss.hDS;
  i_ba.anyelig_idx = p_tx.idx_hm_elig(t);

  internal::run_calculate_incidence_rate<ModelVariant>(
    t, pars, state_curr, state_next, intermediate
  );

  for (int hiv_step = 0; hiv_step < p_op.hts_per_year; ++hiv_step) {
    i_ba.grad.setZero();
    i_ba.gradART.setZero();
    i_ba.p_hiv_deaths_age_sex.setZero();

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
