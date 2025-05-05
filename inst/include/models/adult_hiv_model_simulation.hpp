#pragma once

#include "../options.hpp"
#include "../generated/config_mixer.hpp"

namespace leapfrog {
namespace internal {

template<typename Config>
concept AdultHivModelSimulationEnabled = RunDemographicProjection<Config> && RunHivSimulation<Config>;

template<typename Config>
struct AdultHivModelSimulation {
  AdultHivModelSimulation(...) {};
};

template<AdultHivModelSimulationEnabled Config>
struct AdultHivModelSimulation<Config> {
  using real_type = typename Config::real_type;
  using ModelVariant = typename Config::ModelVariant;
  using SS = Config::SS;
  using Pars = Config::Pars;
  using State = Config::State;
  using Intermediate = Config::Intermediate;
  using Args = Config::Args;

  // private members of this struct
  private:
  // state space
  static constexpr int NS = SS::NS;
  static constexpr int pAG = SS::pAG;
  static constexpr int hDS = SS::hDS;
  static constexpr int hTS = SS::hTS;
  static constexpr int hAG = SS::hAG;
  static constexpr auto hAG_span = SS::hAG_span;
  static constexpr int PROJPERIOD_MIDYEAR = SS::PROJPERIOD_MIDYEAR;
  static constexpr int MALE = SS::MALE;
  static constexpr int FEMALE = SS::FEMALE;
  static constexpr int ART0MOS = SS::ART0MOS;
  static constexpr int p_idx_hiv_first_adult = SS::p_idx_hiv_first_adult;
  static constexpr int adult_incidence_first_age_group = SS::adult_incidence_first_age_group;
  static constexpr int pAG_INCIDPOP = SS::pAG_INCIDPOP;
  static constexpr int hIDX_15PLUS = SS::hIDX_15PLUS;

  // function args
  int t;
  const Pars& pars;
  const State& state_curr;
  State& state_next;
  Intermediate& intermediate;
  const Options<real_type>& opts;

  // only exposing the constructor and some methods
  public:
  AdultHivModelSimulation(Args& args):
    t(args.t),
    pars(args.pars),
    state_curr(args.state_curr),
    state_next(args.state_next),
    intermediate(args.intermediate),
    opts(args.opts)
  {};

  void run_hiv_model_simulation() {
    const auto& p_ha = pars.ha;
    auto& i_ha = intermediate.ha;

    i_ha.everARTelig_idx = p_ha.idx_hm_elig(t) < hDS ? p_ha.idx_hm_elig(t) : hDS;
    i_ha.anyelig_idx = p_ha.idx_hm_elig(t);
    run_calculate_incidence_rate();

    for (int hiv_step = 0; hiv_step < opts.hts_per_year; ++hiv_step) {
      i_ha.grad.setZero();
      i_ha.gradART.setZero();
      i_ha.p_hiv_deaths_age_sex.setZero();

      run_disease_progression_and_mortality(hiv_step);
      run_new_p_infections(hiv_step);
      run_new_hiv_p_infections(hiv_step);

      if (t >= opts.ts_art_start) {
        run_art_progression_and_mortality(hiv_step);
        run_h_art_initiation(hiv_step);
        run_update_art_stratification(hiv_step);
      }

      run_update_hiv_stratification(hiv_step);
      run_remove_p_hiv_deaths(hiv_step);
    }
  };

  // private methods that we don't want people to call
  private:
  void distribute_rate_over_sexes() {
    const auto& p_ha = pars.ha;
    auto& i_ha = intermediate.ha;

    real_type denominator = i_ha.hiv_neg_aggregate(MALE) +
                            p_ha.relative_risk_sex(t) * i_ha.hiv_neg_aggregate(FEMALE);
    real_type total_neg = i_ha.hiv_neg_aggregate(MALE) + i_ha.hiv_neg_aggregate(FEMALE);
    i_ha.rate_sex(MALE) = p_ha.total_rate(t) * total_neg / denominator;
    i_ha.rate_sex(FEMALE) = i_ha.rate_sex(MALE) * p_ha.relative_risk_sex(t);
  };

  void run_calculate_incidence_rate() {
    const auto& c_dp = state_curr.dp;
    const auto& c_ha = state_curr.ha;
    auto& i_ha = intermediate.ha;

    const auto adult_incid_first_age_group = adult_incidence_first_age_group;
    const auto adult_incid_last_age_group = adult_incid_first_age_group + pAG_INCIDPOP;

    for (int g = 0; g < NS; ++g) {
      for (int a = adult_incid_first_age_group; a < adult_incid_last_age_group; ++a) {
        i_ha.hiv_neg_aggregate(g) += c_dp.p_total_pop(a, g) - c_ha.p_hiv_pop(a, g);
      }
    }

    distribute_rate_over_sexes();
  };

  void run_disease_progression_and_mortality(int hiv_step) {
    const auto& p_ha = pars.ha;
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;

    for (int g = 0; g < NS; ++g) {
      for (int ha = 0; ha < hAG; ++ha) {
        for (int hm = 0; hm < hDS; ++hm) {
          i_ha.cd4mx_scale = 1.0;
          if (p_ha.scale_cd4_mortality && t >= opts.ts_art_start &&
              hm >= i_ha.everARTelig_idx && n_ha.h_hiv_adult(hm, ha, g) > 0.0) {
            i_ha.artpop_hahm = 0.0;
            for (int hu = 0; hu < hTS; ++hu) {
              i_ha.artpop_hahm += n_ha.h_art_adult(hu, hm, ha, g);
            }
            i_ha.cd4mx_scale = n_ha.h_hiv_adult(hm, ha, g) /
                              (n_ha.h_hiv_adult(hm, ha, g) + i_ha.artpop_hahm);
          }

          i_ha.deaths = i_ha.cd4mx_scale * p_ha.cd4_mortality(hm, ha, g) * n_ha.h_hiv_adult(hm, ha, g);
          i_ha.p_hiv_deaths_age_sex(ha, g) += opts.dt * i_ha.deaths;
          n_ha.h_hiv_deaths_no_art(hm, ha, g) += opts.dt * i_ha.deaths;
          i_ha.grad(hm, ha, g) = -i_ha.deaths;
        }

        for (int hm = 1; hm < hDS; ++hm) {
          const auto hiv_adults_progressing_cd4_stage = p_ha.cd4_progression(hm - 1, ha, g) * n_ha.h_hiv_adult(hm - 1, ha, g);
          i_ha.grad(hm - 1, ha, g) -= hiv_adults_progressing_cd4_stage;
          i_ha.grad(hm, ha, g) += hiv_adults_progressing_cd4_stage;
        }
      }
    }
  };

  void run_new_p_infections(int hiv_step) {
    const auto& p_ha = pars.ha;
    auto& n_ha = state_next.ha;
    auto& n_dp = state_next.dp;
    auto& i_ha = intermediate.ha;

    const auto adult_incid_first_age_group = adult_incidence_first_age_group;
    const auto adult_incid_last_age_group = adult_incid_first_age_group + pAG_INCIDPOP;

    for (int g = 0; g < NS; ++g) {
      i_ha.hiv_negative_pop.setZero();
      i_ha.Xhivn_incagerr = 0.0;

      for (int a = adult_incid_first_age_group; a < pAG; ++a) {
        i_ha.hiv_negative_pop(a) = n_dp.p_total_pop(a, g) - n_ha.p_hiv_pop(a, g);
      }

      for (int a = adult_incid_first_age_group; a < adult_incid_last_age_group; ++a) {
        i_ha.Xhivn_incagerr += p_ha.relative_risk_age(a - adult_incid_first_age_group, g, t) *
                               i_ha.hiv_negative_pop(a);
      }

      for (int a = adult_incid_first_age_group; a < pAG; ++a) {
        i_ha.p_infections_ts(a, g) = i_ha.hiv_negative_pop(a) *
                                     i_ha.rate_sex(g) *
                                     p_ha.relative_risk_age(a - adult_incid_first_age_group, g, t) *
                                     i_ha.hiv_neg_aggregate(g) /
                                     i_ha.Xhivn_incagerr;
      }
    }
  };

  void run_new_hiv_p_infections(int hiv_step) {
    const auto& p_ha = pars.ha;
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;

    for (int g = 0; g < NS; g++) {
      int a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {
        i_ha.p_infections_ha = 0.0;
        for (int i = 0; i < hAG_span[ha]; i++, a++) {
          i_ha.p_infections_a = i_ha.p_infections_ts(a, g);
          i_ha.p_infections_ha += i_ha.p_infections_a;
          const auto new_infections = opts.dt * i_ha.p_infections_a;
          n_ha.p_infections(a, g) += new_infections;
          n_ha.p_hiv_pop(a, g) += new_infections;
        }

        // add p_infections to grad hivpop
        for (int hm = 0; hm < hDS; ++hm) {
          i_ha.grad(hm, ha, g) += i_ha.p_infections_ha * p_ha.cd4_initial_distribution(hm, ha, g);
        }
      }
    }
  };

  void run_art_progression_and_mortality(int hiv_step) {
    const auto& p_ha = pars.ha;
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;

    for (int g = 0; g < NS; ++g) {
      for (int ha = 0; ha < hAG; ++ha) {
        for (int hm = i_ha.everARTelig_idx; hm < hDS; ++hm) {
          for (int hu = 0; hu < hTS; ++hu) {
            i_ha.deaths_art = p_ha.mortality(hu, hm, ha, g) *
                              p_ha.mortality_time_rate_ratio(hu, t) *
                              n_ha.h_art_adult(hu, hm, ha, g);
            const auto new_hiv_deaths_art = opts.dt * i_ha.deaths_art;
            i_ha.p_hiv_deaths_age_sex(ha, g) += new_hiv_deaths_art;
            n_ha.h_hiv_deaths_art(hu, hm, ha, g) += new_hiv_deaths_art;
            i_ha.gradART(hu, hm, ha, g) = -i_ha.deaths_art;
          }

          for (int hu = 0; hu < (hTS - 1); ++hu) {
            const auto art_adults_progressing_treatment_stage = n_ha.h_art_adult(hu, hm, ha, g) / p_ha.h_art_stage_dur(hu);
            i_ha.gradART(hu, hm, ha, g) -= art_adults_progressing_treatment_stage;
            i_ha.gradART(hu + 1, hm, ha, g) += art_adults_progressing_treatment_stage;
          }

          // ART dropout
          if (p_ha.dropout_rate(t) > 0) {
            for (int hu = 0; hu < hTS; ++hu) {
              const auto art_adult_dropout = p_ha.dropout_rate(t) * n_ha.h_art_adult(hu, hm, ha, g);
              if (p_ha.dropout_recover_cd4 && hu >= 2 && hm >= 1) {
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
  };

  void run_h_art_initiation(int hiv_step) {
    const auto& p_ha = pars.ha;
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;

    for (int g = 0; g < NS; ++g) {
      i_ha.Xart_15plus = 0.0;

      i_ha.artelig_hm.setZero();
      i_ha.Xartelig_15plus = 0.0;

      i_ha.expect_mort_artelig_hm.setZero();
      i_ha.expect_mort_artelig15plus = 0.0;

      for (int ha = hIDX_15PLUS; ha < hAG; ++ha) {
        for (int hm = i_ha.everARTelig_idx; hm < hDS; ++hm) {
          if (hm >= i_ha.anyelig_idx) {
            // TODO: Implement special population ART eligibility
            real_type prop_elig = 1.0;
            real_type tmp_artelig = prop_elig * n_ha.h_hiv_adult(hm, ha, g);
            i_ha.artelig_hahm(hm, ha - hIDX_15PLUS) = tmp_artelig;
            i_ha.artelig_hm(hm) += tmp_artelig;
            i_ha.Xartelig_15plus += tmp_artelig;

            real_type tmp_expect_mort = p_ha.cd4_mortality(hm, ha, g) * i_ha.artelig_hahm(hm, ha - hIDX_15PLUS);
            i_ha.expect_mort_artelig_hm(hm) += tmp_expect_mort;
            i_ha.expect_mort_artelig15plus += tmp_expect_mort;
          }

          for (int hu = 0; hu < hTS; ++hu) {
            i_ha.Xart_15plus += n_ha.h_art_adult(hu, hm, ha, g) +
                                opts.dt * i_ha.gradART(hu, hm, ha, g);
          }
        }
      }

      // calculate number on ART at end of ts, based on number or percent
      real_type art_interp_w = opts.dt * (hiv_step + 1.0);
      if (opts.proj_period_int == PROJPERIOD_MIDYEAR && art_interp_w < 0.5) {
        if (!p_ha.adults_on_art_is_percent(g, t - 2) && !p_ha.adults_on_art_is_percent(g, t - 1)) {
          // case: both values are numbers
          i_ha.artnum_hts = (0.5 - art_interp_w) * p_ha.adults_on_art(g, t - 2) +
                            (art_interp_w + 0.5) * p_ha.adults_on_art(g, t - 1);
        } else if (p_ha.adults_on_art_is_percent(g, t - 2) && p_ha.adults_on_art_is_percent(g, t - 1)) {
          // case: both values are percentages
          i_ha.artcov_hts = (0.5 - art_interp_w) * p_ha.adults_on_art(g, t - 2) +
                            (art_interp_w + 0.5) * p_ha.adults_on_art(g, t - 1);
          i_ha.artnum_hts = i_ha.artcov_hts * (i_ha.Xart_15plus + i_ha.Xartelig_15plus);
        } else if (!p_ha.adults_on_art_is_percent(g, t - 2) && p_ha.adults_on_art_is_percent(g, t - 1)) {
          // case: value is percentage only at time t - 1
          // transition from number to percentage
          i_ha.curr_coverage = i_ha.Xart_15plus / (i_ha.Xart_15plus + i_ha.Xartelig_15plus);
          i_ha.artcov_hts = i_ha.curr_coverage +
                            (p_ha.adults_on_art(g, t - 1) - i_ha.curr_coverage) *
                            opts.dt / (0.5 - opts.dt * hiv_step);
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

        if (opts.proj_period_int == PROJPERIOD_MIDYEAR) {
          art_interp_w -= 0.5;
        }

        if (!p_ha.adults_on_art_is_percent(g, t - 1) && !p_ha.adults_on_art_is_percent(g, t)) {
          // case: both values are numbers
          i_ha.artnum_hts = (1.0 - art_interp_w) * p_ha.adults_on_art(g, t - 1) +
                            art_interp_w * p_ha.adults_on_art(g, t);
        } else if (p_ha.adults_on_art_is_percent(g, t - 1) && p_ha.adults_on_art_is_percent(g, t)) {
          // case: both values are percentages
          i_ha.artcov_hts = (1.0 - art_interp_w) * p_ha.adults_on_art(g, t - 1) +
                            art_interp_w * p_ha.adults_on_art(g, t);
          // transition to number
          i_ha.artnum_hts = i_ha.artcov_hts * (i_ha.Xart_15plus + i_ha.Xartelig_15plus);
        } else if (!p_ha.adults_on_art_is_percent(g, t - 1) && p_ha.adults_on_art_is_percent(g, t)) {
          // case: value is percentage only at time t
          // transition from number to percentage
          i_ha.curr_coverage = i_ha.Xart_15plus / (i_ha.Xart_15plus + i_ha.Xartelig_15plus);
          i_ha.artcov_hts = i_ha.curr_coverage +
                            (p_ha.adults_on_art(g, t) - i_ha.curr_coverage) *
                            opts.dt / (1.0 - art_interp_w + opts.dt);
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
      for (int hm = i_ha.anyelig_idx; hm < hDS; ++hm) {
        auto eligibility_by_stage = (1.0 - p_ha.initiation_mortality_weight) *
                                    i_ha.artelig_hm(hm) /
                                    i_ha.Xartelig_15plus;
        auto expected_mortality_by_stage = p_ha.initiation_mortality_weight *
                                          i_ha.expect_mort_artelig_hm(hm) /
                                          i_ha.expect_mort_artelig15plus;
        i_ha.artinit_hm(hm) = i_ha.artinit_hts * (eligibility_by_stage + expected_mortality_by_stage);
      }

      // Step 2: within CD4 category, allocate ART by age proportional to
      // eligibility
      for (int ha = hIDX_15PLUS; ha < hAG; ++ha) {
        for (int hm = i_ha.anyelig_idx; hm < hDS; ++hm) {
          if (i_ha.artelig_hm(hm) > 0.0) {
            i_ha.artinit_hahm = i_ha.artinit_hm(hm) *
                                i_ha.artelig_hahm(hm, ha - hIDX_15PLUS) /
                                i_ha.artelig_hm(hm);
            i_ha.artinit_hahm = std::min(i_ha.artinit_hahm, i_ha.artelig_hahm(hm, ha - hIDX_15PLUS));
            i_ha.artinit_hahm = std::min(i_ha.artinit_hahm,
                                        n_ha.h_hiv_adult(hm, ha, g) + opts.dt * i_ha.grad(hm, ha, g));
            i_ha.grad(hm, ha, g) -= i_ha.artinit_hahm / opts.dt;
            i_ha.gradART(ART0MOS, hm, ha, g) += i_ha.artinit_hahm / opts.dt;
            n_ha.h_art_initiation(hm, ha, g) += i_ha.artinit_hahm;
          }
        }
      }
    }
  };

  void run_update_art_stratification(int hiv_step) {
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;

    for (int g = 0; g < NS; ++g) {
      for (int ha = 0; ha < hAG; ++ha) {
        for (int hm = i_ha.everARTelig_idx; hm < hDS; ++hm) {
          for (int hu = 0; hu < hTS; ++hu) {
            n_ha.h_art_adult(hu, hm, ha, g) += opts.dt * i_ha.gradART(hu, hm, ha, g);
          }
        }
      }
    }
  };

  void run_update_hiv_stratification(int hiv_step) {
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;

    for (int g = 0; g < NS; ++g) {
      for (int ha = 0; ha < hAG; ++ha) {
        for (int hm = 0; hm < hDS; ++hm) {
          n_ha.h_hiv_adult(hm, ha, g) += opts.dt * i_ha.grad(hm, ha, g);
        }
      }
    }
  };

  void run_remove_p_hiv_deaths(int hiv_step) {
    auto& n_dp = state_next.dp;
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;

    for (int g = 0; g < NS; ++g) {
      // sum HIV+ population size in each hivpop age group
      int a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {
        i_ha.hivpop_ha(ha) = 0.0;
        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
          i_ha.hivpop_ha(ha) += n_ha.p_hiv_pop(a, g);
        }
      }

      // remove hivdeaths proportionally to age-distribution within each age group
      a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {
        if (i_ha.hivpop_ha(ha) > 0) {
          i_ha.hivqx_ha = i_ha.p_hiv_deaths_age_sex(ha, g) / i_ha.hivpop_ha(ha);
          for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
            i_ha.hivdeaths_a = n_ha.p_hiv_pop(a, g) * i_ha.hivqx_ha;
            n_ha.p_hiv_deaths(a, g) += i_ha.hivdeaths_a;
            n_dp.p_total_pop(a, g) -= i_ha.hivdeaths_a;
            n_ha.p_hiv_pop(a, g) -= i_ha.hivdeaths_a;
          }
        } else {
          a += hAG_span[ha];
        }
      }
    }
  };
};

}
}
