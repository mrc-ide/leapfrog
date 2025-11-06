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

  void run_hivpop_demographic_projection() {

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

  void run_hivpop_end_year_migration() {
    auto& n_ha = state_next.ha;
    auto& i_dp = intermediate.dp;

    // remove net migration from hiv stratified population
    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < pAG; ++a) {
        n_ha.p_net_migration_hivpop(a, s) = n_ha.p_hivpop(a, s) * i_dp.migration_rate(a, s);
        n_ha.p_hivpop(a, s) += n_ha.p_net_migration_hivpop(a, s);
      }
    }

    // remove net migration from adult stratified population
    for (int s = 0; s < NS; ++s) {
      int a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {
        real_type migration_num_ha = 0.0;
        real_type hivpop_ha_postmig = 0.0;
        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
          hivpop_ha_postmig += n_ha.p_hivpop(a, s);
          migration_num_ha += n_ha.p_net_migration_hivpop(a, s);
        }

        real_type migration_rate = 0.0;
        if (hivpop_ha_postmig > 0.0) {
          migration_rate = migration_num_ha / (hivpop_ha_postmig - migration_num_ha);
        }

        for (int hm = 0; hm < hDS; ++hm) {
          n_ha.h_hivpop(hm, ha, s) *= 1.0 + migration_rate;
          if (t >= opts.ts_art_start) {
            for (int hu = 0; hu < hTS; ++hu) {
              n_ha.h_artpop(hu, hm, ha, s) *= 1.0 + migration_rate;
            }
          }
        }
      }
    }
  };

  void run_hc_hivpop_end_year_migration() {
    static_assert(ModelVariant::run_child_model,
                  "run_hc_hivpop_end_year_migration can only be called for model variants where run_child_model is true");
    static constexpr int hc2_agestart = SS::hc2_agestart;
    static constexpr int hcAG_end = SS::hcAG_end;
    static constexpr int hc1DS = SS::hc1DS;
    static constexpr int hc2DS = SS::hc2DS;
    static constexpr int hTS = SS::hTS;
    static constexpr int hcTT = SS::hcTT;

    auto& n_ha = state_next.ha;
    auto& n_hc = state_next.hc;
    const auto& p_hc = pars.hc;

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < hcAG_end; ++a) {
        real_type hc_migration_num = 0.0;
        real_type hc_hivpop_postmig = n_ha.p_hivpop(a, s);
        hc_migration_num = n_ha.p_net_migration_hivpop(a, s);

        real_type migration_rate = 0.0;
        if (hc_hivpop_postmig > 0.0) {
          migration_rate = hc_migration_num / (hc_hivpop_postmig - hc_migration_num);
        }

        if (a < hc2_agestart) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            for (int cat = 0; cat < hcTT; ++cat) {
              n_hc.hc1_hivpop(hd, cat, a, s) *= 1.0 + migration_rate;
            }
            if (t >= p_hc.hc_art_start) {
              for (int dur = 0; dur < hTS; ++dur) {
                n_hc.hc1_artpop(dur, hd, a, s) *= 1.0 + migration_rate;
              }
            }
          }
        } else {
          for (int hd = 0; hd < hc2DS; ++hd) {
            for (int cat = 0; cat < hcTT; ++cat) {
              n_hc.hc2_hivpop(hd, cat, a - hc2_agestart, s) *= 1.0 + migration_rate;
            }
            if (t >= p_hc.hc_art_start) {
              for (int dur = 0; dur < hTS; ++dur) {
                n_hc.hc2_artpop(dur, hd, a - hc2_agestart, s) *= 1.0 + migration_rate;
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
    for (int s = 0; s < NS; ++s) {
      for (int a = 1; a < pAG; ++a) {
        n_ha.p_deaths_background_hivpop(a, s) = c_ha.p_hivpop(a - 1, s) * (1.0 - p_dp.survival_probability(a, s, t));
        n_ha.p_hivpop(a, s) = c_ha.p_hivpop(a - 1, s);
      }

      // open age group
      n_ha.p_deaths_background_hivpop(pAG - 1, s) += c_ha.p_hivpop(pAG - 1, s) *
                                                   (1.0 - p_dp.survival_probability(pAG, s, t));
      n_ha.p_hivpop(pAG - 1, s) += c_ha.p_hivpop(pAG - 1, s);
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

    for (int s = 0; s < NS; ++s) {
      for (int hm = 0; hm < hDS; ++hm) {
	i_hc.age15_hivpop(hm, s) = 0.0;
	for (int hm_adol = 0; hm_adol < hc2DS; ++hm_adol) {
	  auto age15_hivpop_hm_adol = 0.0;
	  for (int htm = 0; htm < hcTT; ++htm) {
	    age15_hivpop_hm_adol += c_hc.hc2_hivpop(hm_adol, htm, (hc2AG - 1), s);
	  }
	  i_hc.age15_hivpop(hm, s) += age15_hivpop_hm_adol * SS::hc2_to_ha_cd4_transition[hm][hm_adol];
	}

	if (t > p_hc.hc_art_start) {
	  for (int hu = 0; hu < hTS; ++hu) {
	    i_hc.age15_artpop(hu, hm, s) = 0.0;
	    for (int hm_adol = 0; hm_adol < hc2DS; ++hm_adol) {
	      i_hc.age15_artpop(hu, hm, s) += c_hc.hc2_artpop(hu, hm_adol, (hc2AG - 1), s)  * SS::hc2_to_ha_cd4_transition[hm][hm_adol];
	    }
	  }
	}
      } // hm
    } // s
  };
  
  void run_hiv_and_art_stratified_ageing() {
    const auto& c_ha = state_curr.ha;
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;

    // age coarse stratified HIV population
    for (int s = 0; s < NS; ++s) {
      int a = p_idx_hiv_first_adult;
      // Note: loop stops at hAG-1; no one ages out of the open-ended
      // age group
      for (int ha = 0; ha < (hAG - 1); ++ha) {
        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
          i_ha.hiv_age_up_prob(ha, s) += c_ha.p_hivpop(a, s);
        }

        if (i_ha.hiv_age_up_prob(ha, s) > 0) {
          i_ha.hiv_age_up_prob(ha, s) = c_ha.p_hivpop(a - 1, s) / i_ha.hiv_age_up_prob(ha, s);
        } else {
          i_ha.hiv_age_up_prob(ha, s) = 0.0;
        }
      }
    }

    for (int s = 0; s < NS; ++s) {
      for (int ha = 1; ha < hAG; ++ha) {
        for (int hm = 0; hm < hDS; ++hm) {
          n_ha.h_hivpop(hm, ha, s) = (1.0 - i_ha.hiv_age_up_prob(ha, s)) * c_ha.h_hivpop(hm, ha, s);  // age-out
	  n_ha.h_hivpop(hm, ha, s) += i_ha.hiv_age_up_prob(ha - 1, s) * c_ha.h_hivpop(hm, ha - 1, s); // age-in

          if (t > opts.ts_art_start) {
            for (int hu = 0; hu < hTS; ++hu) {
              n_ha.h_artpop(hu, hm, ha, s) = (1.0 - i_ha.hiv_age_up_prob(ha, s)) * c_ha.h_artpop(hu, hm, ha, s);  // age-out
	      n_ha.h_artpop(hu, hm, ha, s) += i_ha.hiv_age_up_prob(ha - 1, s) * c_ha.h_artpop(hu, hm, ha - 1, s); // age-in
            }
	  }
        }
      }

      // age out of ha = 0 group
      int ha = 0;
      for (int hm = 0; hm < hDS; ++hm) {
	n_ha.h_hivpop(hm, ha, s) = (1.0 - i_ha.hiv_age_up_prob(ha, s)) * c_ha.h_hivpop(hm, ha, s);  // age-out
	if (t > opts.ts_art_start) {
	  for (int hu = 0; hu < hTS; ++hu) {
	    n_ha.h_artpop(hu, hm, ha, s) = (1.0 - i_ha.hiv_age_up_prob(ha, s)) * c_ha.h_artpop(hu, hm, ha, s);  // age-out
	  }
	}
      }
      
    // Entrants ageing into adult HIV population
    if constexpr (ModelVariant::run_child_model) {
      auto& i_hc = intermediate.hc;
      const auto& p_hc = pars.hc;

        for (int hm = 0; hm < hDS; ++hm) {
	  n_ha.h_hivpop(hm, 0, s) += i_hc.age15_hivpop(hm, s);

	  if (t > p_hc.hc_art_start) {
	    for (int hu = 0; hu < hTS; ++hu) {
	      if (t > opts.ts_art_start) {
		n_ha.h_artpop(hu, hm, 0, s) += i_hc.age15_artpop(hu, hm, s);
	      } else {
		// If child ART has started, but not yet adult ART has started, put
		// children on ART into the untreated adult HIV population
		n_ha.h_hivpop(hm, 0, s) += i_hc.age15_artpop(hu, hm, s);
	      }
	    }
	  }
	}
      }
    } // s
    // TODO: implement static entrants to adult HIV population here for case when child model not simulated
    
  };

  void run_hiv_and_art_stratified_deaths_and_migration() {
    auto& n_ha = state_next.ha;
    auto& i_ha = intermediate.ha;
    auto& i_dp = intermediate.dp;

    for (int s = 0; s < NS; ++s) {
      int a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {
        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
          i_ha.p_coarse_ages_hivpop(ha, s) += n_ha.p_hivpop(a, s);
        }
      }
    }

    // remove non-HIV deaths and net migration from hiv stratified population
    for (int s = 0; s < NS; ++s) {
      for (int a = 1; a < pAG; ++a) {
        n_ha.p_hivpop(a, s) -= n_ha.p_deaths_background_hivpop(a, s);
        if (opts.proj_period_int == PROJPERIOD_MIDYEAR) {
          n_ha.p_net_migration_hivpop(a, s) = n_ha.p_hivpop(a, s) * i_dp.migration_rate(a, s);
          n_ha.p_hivpop(a, s) += n_ha.p_net_migration_hivpop(a, s);
        }
      }
    }
    
    // remove non-HIV deaths and net migration from adult stratified population
    for (int s = 0; s < NS; ++s) {
      int a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {
        real_type deaths_migrate = 0.0;
        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
          deaths_migrate -= n_ha.p_deaths_background_hivpop(a, s);
          if (opts.proj_period_int == PROJPERIOD_MIDYEAR) {
            deaths_migrate += n_ha.p_net_migration_hivpop(a, s);
          }
        }

        real_type deaths_migrate_rate = 0.0;
        if (i_ha.p_coarse_ages_hivpop(ha, s) > 0) {
          deaths_migrate_rate = deaths_migrate / i_ha.p_coarse_ages_hivpop(ha, s);
        }

        for (int hm = 0; hm < hDS; ++hm) {
          n_ha.h_hivpop(hm, ha, s) *= 1.0 + deaths_migrate_rate;
          if (t > opts.ts_art_start) {
            for (int hu = 0; hu < hTS; ++hu) {
              n_ha.h_artpop(hu, hm, ha, s) *= 1.0 + deaths_migrate_rate;
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
    const auto& p_hc = pars.hc;

    real_type deaths_migrate = 0.0;
    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < hcAG_end; ++a) {
        deaths_migrate -= n_ha.p_deaths_background_hivpop(a, s);
        if (opts.proj_period_int == PROJPERIOD_MIDYEAR) {
          deaths_migrate += n_ha.p_net_migration_hivpop(a, s);
        }
        if(n_ha.p_hivpop(a, s) > 0){
          deaths_migrate /= n_ha.p_hivpop(a, s);
        }
        if (a < hc2_agestart) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            for (int cat = 0; cat < hcTT; ++cat) {
              n_hc.hc1_hivpop(hd, cat, a, s) *= 1.0 + deaths_migrate;
            }
            if (t > p_hc.hc_art_start) {
              for (int dur = 0; dur < hTS; ++dur) {
                n_hc.hc1_artpop(dur, hd, a, s) *= 1.0 + deaths_migrate;
              }
            }
          }
          } else {
            for (int hd = 0; hd < hc2DS; ++hd) {
              for (int cat = 0; cat < hcTT; ++cat) {
                n_hc.hc2_hivpop(hd, cat, a - hc2_agestart, s) *= 1.0 + deaths_migrate;
              }
              if (t > p_hc.hc_art_start) {
                for (int dur = 0; dur < hTS; ++dur) {
                  n_hc.hc2_artpop(dur, hd, a - hc2_agestart, s) *= 1.0 + deaths_migrate;
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
