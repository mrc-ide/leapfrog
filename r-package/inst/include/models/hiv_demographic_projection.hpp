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
    auto& n_ha = state_next.ha;
    auto& c_ha = state_curr.ha;

    run_hiv_ageing_and_mortality();
    if constexpr (ModelVariant::run_child_model) {
      run_age_15_entrants();
    }
    run_hiv_and_art_stratified_ageing();
    run_hiv_and_art_stratified_deaths_and_migration();
    if constexpr (ModelVariant::run_child_model) {
      run_hc_hiv_and_art_stratified_deaths_and_migration();
    }

  };

  void run_hiv_pop_end_year_migration() {
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;
    auto& i_dp = intermediate.dp;

    // remove net migration from hiv stratified population
    for (int g = 0; g < NS; ++g) {
      for (int a = 0; a < pAG; ++a) {
        n_ha.p_net_migration_hivpop(a, g) = n_ha.p_hiv_pop(a, g) * i_dp.migration_rate(a, g);
        n_ha.p_hiv_pop(a, g) += n_ha.p_net_migration_hivpop(a, g);
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
          migration_num_ha += n_ha.p_net_migration_hivpop(a, g);
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

  void run_hc_hiv_pop_end_year_migration() {
    static_assert(ModelVariant::run_child_model,
                  "run_age_15_entrants can only be called for model variants where run_child_model is true");
    static constexpr int hc2_agestart = SS::hc2_agestart;
    static constexpr int hcAG_end = SS::hcAG_end;
    static constexpr int hc1DS = SS::hc1DS;
    static constexpr int hc2DS = SS::hc2DS;
    static constexpr int hTS = SS::hTS;
    static constexpr int hcTT = SS::hcTT;

    auto& n_ha = state_next.ha;
    auto& n_hc = state_next.hc;
    auto& i_ha = intermediate.ha;
    auto& i_dp = intermediate.dp;
    const auto& p_hc = pars.hc;

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < hcAG_end; ++a) {
        real_type hc_migration_num = 0.0;
        real_type hc_hivpop_postmig = n_ha.p_hiv_pop(a, s);
        hc_migration_num = n_ha.p_net_migration_hivpop(a, s);

        real_type migration_rate = 0.0;
        if (hc_hivpop_postmig > 0.0) {
          migration_rate = hc_migration_num / (hc_hivpop_postmig - hc_migration_num);
        }

        if (a < hc2_agestart) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            for (int cat = 0; cat < hcTT; ++cat) {
              n_hc.hc1_hiv_pop(hd, cat, a, s) *= 1.0 + migration_rate;
            }
            if (t >= p_hc.hc_art_start) {
              for (int dur = 0; dur < hTS; ++dur) {
                n_hc.hc1_art_pop(dur, hd, a, s) *= 1.0 + migration_rate;
              }
            }
          }
        } else {
          for (int hd = 0; hd < hc2DS; ++hd) {
            for (int cat = 0; cat < hcTT; ++cat) {
              n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) *= 1.0 + migration_rate;
            }
            if (t >= p_hc.hc_art_start) {
              for (int dur = 0; dur < hTS; ++dur) {
                n_hc.hc2_art_pop(dur, hd, a - hc2_agestart, s) *= 1.0 + migration_rate;
              }
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

    const auto& p_hc = pars.hc;
    const auto& c_hc = state_curr.hc;
    auto& i_hc = intermediate.hc;

    for (int g = 0; g < NS; ++g) {
      for (int hm = 0; hm < hDS; ++hm) {
	i_hc.age15_hiv_pop(hm, g) = 0.0;
	for (int hm_adol = 0; hm_adol < hc2DS; ++hm_adol) {
	  auto age15_hivpop_hm_adol = 0.0;
	  for (int htm = 0; htm < hcTT; ++htm) {
	    age15_hivpop_hm_adol += c_hc.hc2_hiv_pop(hm_adol, htm, (hc2AG - 1), g);
	  }
	  i_hc.age15_hiv_pop(hm, g) += age15_hivpop_hm_adol * SS::hc2_to_ha_cd4_transition[hm][hm_adol];
	}

	if (t > p_hc.hc_art_start) {
	  for (int hu = 0; hu < hTS; ++hu) {
	    i_hc.age15_art_pop(hu, hm, g) = 0.0;
	      for (int hm_adol = 0; hm_adol < hc2DS; ++hm_adol) {
		i_hc.age15_art_pop(hu, hm, g) += c_hc.hc2_art_pop(hu, hm_adol, (hc2AG - 1), g)  * SS::hc2_to_ha_cd4_transition[hm][hm_adol];
	      }
	  }
	}
	
      } // hm
    } // g
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
          n_ha.h_hiv_adult(hm, ha, g) = (1.0 - i_ha.hiv_age_up_prob(ha, g)) * c_ha.h_hiv_adult(hm, ha, g);  // age-out
	  n_ha.h_hiv_adult(hm, ha, g) += i_ha.hiv_age_up_prob(ha - 1, g) * c_ha.h_hiv_adult(hm, ha - 1, g); // age-in
          if (t > opts.ts_art_start) {
            for (int hu = 0; hu < hTS; ++hu) {
              n_ha.h_art_adult(hu, hm, ha, g) = (1.0 - i_ha.hiv_age_up_prob(ha, g)) * c_ha.h_art_adult(hu, hm, ha, g);  // age-out
	      n_ha.h_art_adult(hu, hm, ha, g) += i_ha.hiv_age_up_prob(ha - 1, g) * c_ha.h_art_adult(hu, hm, ha - 1, g); // age-in
            }
	  }
        }
      }

      // age out of ha = 0 group
      int ha = 0;
      for (int hm = 0; hm < hDS; ++hm) {
	n_ha.h_hiv_adult(hm, ha, g) = (1.0 - i_ha.hiv_age_up_prob(ha, g)) * c_ha.h_hiv_adult(hm, ha, g);  // age-out
	if (t > opts.ts_art_start) {
	  for (int hu = 0; hu < hTS; ++hu) {
	    n_ha.h_art_adult(hu, hm, ha, g) = (1.0 - i_ha.hiv_age_up_prob(ha, g)) * c_ha.h_art_adult(hu, hm, ha, g);  // age-out
	  }
	}
      }
      
    // Entrants ageing into adult HIV population
    if constexpr (ModelVariant::run_child_model) {
      auto& i_hc = intermediate.hc;
      const auto& p_hc = pars.hc;

        for (int hm = 0; hm < hDS; ++hm) {
	  n_ha.h_hiv_adult(hm, 0, g) += i_hc.age15_hiv_pop(hm, g);

	  if (t > p_hc.hc_art_start) {
	    for (int hu = 0; hu < hTS; ++hu) {
	      if (t > opts.ts_art_start) {
		n_ha.h_art_adult(hu, hm, 0, g) += i_hc.age15_art_pop(hu, hm, g);
	      } else {
		// If child ART has started, but not yet adult ART has started, put
		// children on ART into the untreated adult HIV population
		n_ha.h_hiv_adult(hm, 0, g) += i_hc.age15_art_pop(hu, hm, g);
	      }
	    }
	  }
	}
      }
    }
    // TODO: implement static entrants to adult HIV population here for case when child model not simulated
    
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
          n_ha.p_net_migration_hivpop(a, g) = n_ha.p_hiv_pop(a, g) * i_dp.migration_rate(a, g);
          n_ha.p_hiv_pop(a, g) += n_ha.p_net_migration_hivpop(a, g);
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
            deaths_migrate += n_ha.p_net_migration_hivpop(a, g);
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

  void run_hc_hiv_and_art_stratified_deaths_and_migration() {
    static constexpr int hc2_agestart = SS::hc2_agestart;
    static constexpr int hcAG_end = SS::hcAG_end;
    static constexpr int hc1DS = SS::hc1DS;
    static constexpr int hc2DS = SS::hc2DS;
    static constexpr int hTS = SS::hTS;
    static constexpr int hcTT = SS::hcTT;

    auto& n_ha = state_next.ha;
    auto& n_hc = state_next.hc;
    auto& i_ha = intermediate.ha;
    auto& i_dp = intermediate.dp;
    const auto& p_hc = pars.hc;

    real_type deaths_migrate = 0.0;
    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < hcAG_end; ++a) {
        deaths_migrate -= n_ha.p_hiv_pop_background_deaths(a, s);
        if (opts.proj_period_int == PROJPERIOD_MIDYEAR) {
          deaths_migrate += n_ha.p_net_migration_hivpop(a, s);
        }
        if(n_ha.p_hiv_pop(a, s) > 0){
          deaths_migrate /= n_ha.p_hiv_pop(a, s);
        }
        if (a < hc2_agestart) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            for (int cat = 0; cat < hcTT; ++cat) {
              n_hc.hc1_hiv_pop(hd, cat, a, s) *= 1.0 + deaths_migrate;
            }
            if (t > p_hc.hc_art_start) {
              for (int dur = 0; dur < hTS; ++dur) {
                n_hc.hc1_art_pop(dur, hd, a, s) *= 1.0 + deaths_migrate;
              }
            }
          }
          } else {
            for (int hd = 0; hd < hc2DS; ++hd) {
              for (int cat = 0; cat < hcTT; ++cat) {
                n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) *= 1.0 + deaths_migrate;
              }
              if (t > p_hc.hc_art_start) {
                for (int dur = 0; dur < hTS; ++dur) {
                  n_hc.hc2_art_pop(dur, hd, a - hc2_agestart, s) *= 1.0 + deaths_migrate;
                }
              }
            }
          }

      }
    }
  };

};
}
}
