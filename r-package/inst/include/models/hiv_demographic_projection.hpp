#pragma once

#include "../options.hpp"
#include "../generated/config_mixer.hpp"

namespace leapfrog {
namespace internal {

template<typename Config>
concept HivDemographicProjectionEnabled = RunDemographicProjection<Config> && RunHivSimulation<Config>;

template<typename Config>
struct HivDemographicProjection {
  HivDemographicProjection(...) {};
};

template<HivDemographicProjectionEnabled Config>
struct HivDemographicProjection<Config> {
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
  static constexpr int p_idx_hiv_first_adult = SS::p_idx_hiv_first_adult;

  // function args
  int t;
  const Pars& pars;
  const State& state_curr;
  State& state_next;
  Intermediate& intermediate;
  const Options<real_type>& opts;

  // only exposing the constructor and some methods
  public:
  HivDemographicProjection(Args& args):
    t(args.t),
    pars(args.pars),
    state_curr(args.state_curr),
    state_next(args.state_next),
    intermediate(args.intermediate),
    opts(args.opts)
  {};

  void run_hiv_pop_demographic_projection() {
    run_hiv_ageing_and_mortality();

    if constexpr (ModelVariant::run_child_model) {
      run_age_15_entrants();
    }

    run_hiv_and_art_stratified_ageing();
    run_hiv_and_art_stratified_deaths_and_migration();
  };

  void run_hiv_pop_end_year_migration() {
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;
    auto& i_dp = intermediate.dp;

    // remove net migration from hiv stratified population
    for (int g = 0; g < NS; ++g) {
      for (int a = 0; a < pAG; ++a) {
        i_ha.hiv_net_migration(a, g) = n_ha.p_hiv_pop(a, g) * i_dp.migration_rate(a, g);
        n_ha.p_hiv_pop(a, g) += i_ha.hiv_net_migration(a, g);
      }
    }

    // remove net migration from adult stratified population
    for (int g = 0; g < NS; ++g) {
      int a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {
        real_type migration_num_ha = 0.0;
        real_type hivpop_ha_postmig = 0.0;
        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
          hivpop_ha_postmig += n_ha.p_hiv_pop(a, g);
          migration_num_ha += i_ha.hiv_net_migration(a, g);
        }

        real_type migration_rate = 0.0;
        if (hivpop_ha_postmig > 0.0) {
          migration_rate = migration_num_ha / (hivpop_ha_postmig - migration_num_ha);
        }

        for (int hm = 0; hm < hDS; ++hm) {
          n_ha.h_hiv_adult(hm, ha, g) *= 1.0 + migration_rate;
          if (t >= opts.ts_art_start) {
            for (int hu = 0; hu < hTS; ++hu) {
              n_ha.h_art_adult(hu, hm, ha, g) *= 1.0 + migration_rate;
            }
          }
        }
      }
    }
  };

  // private methods that we don't want people to call
  private:
  void run_hiv_ageing_and_mortality() {
    const auto& p_dp = pars.dp;
    const auto& c_ha = state_curr.ha;
    auto& n_ha = state_next.ha;

    // Non-hiv deaths
    for (int g = 0; g < NS; ++g) {
      for (int a = 1; a < pAG; ++a) {
        n_ha.p_hiv_pop_background_deaths(a, g) = c_ha.p_hiv_pop(a - 1, g) * (1.0 - p_dp.survival_probability(a, g, t));
        n_ha.p_hiv_pop(a, g) = c_ha.p_hiv_pop(a - 1, g);
      }

      // open age group
      n_ha.p_hiv_pop_background_deaths(pAG - 1, g) += c_ha.p_hiv_pop(pAG - 1, g) *
                                                   (1.0 - p_dp.survival_probability(pAG, g, t));
      n_ha.p_hiv_pop(pAG - 1, g) += c_ha.p_hiv_pop(pAG - 1, g);
    }
  };

  void run_age_15_entrants() {
    static_assert(ModelVariant::run_child_model,
                  "run_age_15_entrants can only be called for model variants where run_child_model is true");
    constexpr int hcTT = SS::hcTT;
    constexpr int hc2AG = SS::hc2AG;
    constexpr int hc2DS = SS::hc2DS;

    const auto& c_hc = state_curr.hc;
    auto& i_hc = intermediate.hc;

    for (int g = 0; g < NS; ++g) {
      for (int hm = 0; hm < hc2DS; ++hm) {
        for (int htm = 0; htm < hcTT; ++htm) {
          i_hc.age15_hiv_pop(hm, g) += c_hc.hc2_hiv_pop(hm, htm, (hc2AG - 1), g);
        }
      }
    }
    for (int g = 0; g < NS; ++g) {
      for (int hm = 0; hm < hc2DS; ++hm) {
        for (int hu = 0; hu < hTS; ++hu) {
          i_hc.age15_art_pop(hu, hm, g) += c_hc.hc2_art_pop(hu, hm, (hc2AG - 1), g);
        }
      }
    }
  };

  void run_hiv_and_art_stratified_ageing() {
    const auto& c_ha = state_curr.ha;
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;

    // age coarse stratified HIV population
    for (int g = 0; g < NS; ++g) {
      int a = p_idx_hiv_first_adult;
      // Note: loop stops at hAG-1; no one ages out of the open-ended
      // age group
      for (int ha = 0; ha < (hAG - 1); ++ha) {
        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
          i_ha.hiv_age_up_prob(ha, g) += c_ha.p_hiv_pop(a, g);
        }

        if (i_ha.hiv_age_up_prob(ha, g) > 0) {
          i_ha.hiv_age_up_prob(ha, g) = c_ha.p_hiv_pop(a - 1, g) / i_ha.hiv_age_up_prob(ha, g);
        } else {
          i_ha.hiv_age_up_prob(ha, g) = 0.0;
        }
      }
    }


    for (int g = 0; g < NS; ++g) {
      for (int ha = 1; ha < hAG; ++ha) {
        for (int hm = 0; hm < hDS; ++hm) {
          n_ha.h_hiv_adult(hm, ha, g) = ((1.0 - i_ha.hiv_age_up_prob(ha, g)) * c_ha.h_hiv_adult(hm, ha, g)) +
                                        (i_ha.hiv_age_up_prob(ha - 1, g) * c_ha.h_hiv_adult(hm, ha - 1, g));
          if (t > opts.ts_art_start)
            for (int hu = 0; hu < hTS; ++hu) {
              n_ha.h_art_adult(hu, hm, ha, g) = ((1.0 - i_ha.hiv_age_up_prob(ha, g)) * c_ha.h_art_adult(hu, hm, ha, g)) +
                                                (i_ha.hiv_age_up_prob(ha - 1, g) * c_ha.h_art_adult(hu, hm, ha - 1, g));
            }
        }
      }
    }

    if constexpr (ModelVariant::run_child_model) {
      constexpr int hc2DS = SS::hc2DS;
      constexpr auto adult_cd4_dist = SS::adult_cd4_dist;

      auto& i_hc = intermediate.hc;

      for (int g = 0; g < NS; ++g) {
        for (int hm = 0; hm < hDS; ++hm) {
          for (int hm_adol = 0; hm_adol < hc2DS; ++hm_adol){
            n_ha.h_hiv_adult(hm, 0, g) += i_hc.age15_hiv_pop(hm_adol, g) * adult_cd4_dist[hm][hm_adol];
            if ((t > opts.ts_art_start)) {
              for (int hu = 0; hu < hTS; ++hu) {
                n_ha.h_art_adult(hu,hm, 0, g) += i_hc.age15_art_pop(hu, hm_adol, g) * adult_cd4_dist[hm][hm_adol];
              }
            }
          }
        }
      }
    } else {
      for (int g = 0; g < NS; ++g) {
        for (int hm = 0; hm < hDS; ++hm) {
          n_ha.h_hiv_adult(hm, 0, g) = (1.0 - i_ha.hiv_age_up_prob(0, g)) * c_ha.h_hiv_adult(hm, 0, g);
          if (t > opts.ts_art_start) {
            for (int hu = 0; hu < hTS; ++hu) {
              n_ha.h_art_adult(hu, hm, 0, g) = (1.0 - i_ha.hiv_age_up_prob(0, g)) * c_ha.h_art_adult(hu, hm, 0, g);
            }
          }
        }
      }
    }
  };

  void run_hiv_and_art_stratified_deaths_and_migration() {
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;
    auto& i_dp = intermediate.dp;

    for (int g = 0; g < NS; ++g) {
      int a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {
        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
          i_ha.p_hiv_pop_coarse_ages(ha, g) += n_ha.p_hiv_pop(a, g);
        }
      }
    }

    // remove non-HIV deaths and net migration from hiv stratified population
    for (int g = 0; g < NS; ++g) {
      for (int a = 1; a < pAG; ++a) {
        n_ha.p_hiv_pop(a, g) -= n_ha.p_hiv_pop_background_deaths(a, g);
        if (opts.proj_period_int == PROJPERIOD_MIDYEAR) {
          i_ha.hiv_net_migration(a, g) = n_ha.p_hiv_pop(a, g) * i_dp.migration_rate(a, g);
          n_ha.p_hiv_pop(a, g) += i_ha.hiv_net_migration(a, g);
        }
      }
    }

    // remove non-HIV deaths and net migration from adult stratified population
    for (int g = 0; g < NS; ++g) {
      int a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {
        real_type deaths_migrate = 0.0;
        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
          deaths_migrate -= n_ha.p_hiv_pop_background_deaths(a, g);
          if (opts.proj_period_int == PROJPERIOD_MIDYEAR) {
            deaths_migrate += i_ha.hiv_net_migration(a, g);
          }
        }

        real_type deaths_migrate_rate = 0.0;
        if (i_ha.p_hiv_pop_coarse_ages(ha, g) > 0) {
          deaths_migrate_rate = deaths_migrate / i_ha.p_hiv_pop_coarse_ages(ha, g);
        }

        for (int hm = 0; hm < hDS; ++hm) {
          n_ha.h_hiv_adult(hm, ha, g) *= 1.0 + deaths_migrate_rate;
          if (t > opts.ts_art_start) {
            for (int hu = 0; hu < hTS; ++hu) {
              n_ha.h_art_adult(hu, hm, ha, g) *= 1.0 + deaths_migrate_rate;
            }
          }
        }
      }
    }
  };
};

}
}
