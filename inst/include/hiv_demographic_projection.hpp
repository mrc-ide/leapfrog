#pragma once

#include <iostream>
#include "intermediate_data.hpp"
#include "generated/state_types.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void run_hiv_ageing_and_mortality(int t,
                                  const Parameters<ModelVariant, real_type> &pars,
                                  const State<ModelVariant, real_type> &state_curr,
                                  State<ModelVariant, real_type> &state_next,
                                  IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_dm = pars.dp.demography;
  const auto& c_ha = state_curr.hiv;
  auto& n_ha = state_next.hiv;

  // Non-hiv deaths
  for (int g = 0; g < ss_d.NS; ++g) {
    for (int a = 1; a < ss_d.pAG; ++a) {
      n_ha.p_hiv_pop_natural_deaths(a, g) = c_ha.p_hiv_pop(a - 1, g) * (1.0 - p_dm.survival_probability(a, g, t));
      n_ha.p_hiv_pop(a, g) = c_ha.p_hiv_pop(a - 1, g);
    }

    // open age group
    n_ha.p_hiv_pop_natural_deaths(ss_d.pAG - 1, g) += c_ha.p_hiv_pop(ss_d.pAG - 1, g) *
                                                      (1.0 - p_dm.survival_probability(ss_d.pAG, g, t));
    n_ha.p_hiv_pop(ss_d.pAG - 1, g) += c_ha.p_hiv_pop(ss_d.pAG - 1, g);
  }
}

template<typename ModelVariant, typename real_type>
void run_age_15_entrants(int t,
                         const Parameters<ModelVariant, real_type> &pars,
                         const State<ModelVariant, real_type> &state_curr,
                         State<ModelVariant, real_type> &state_next,
                         IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>();
  constexpr auto ss_h = ss.hiv;
  constexpr auto ss_d = ss.dp;
  constexpr auto ss_c = ss.children;
  const auto& c_hc = state_curr.children;
  auto& i_hc = intermediate.children;

  for (int g = 0; g < ss_d.NS; ++g) {
    for (int hm = 0; hm < ss_h.hDS; ++hm) {
      for (int htm = 0; htm < ss_c.hcTT; ++htm) {
        i_hc.age15_hiv_pop(hm, g) += c_hc.hc2_hiv_pop(hm, htm, (ss_c.hc2AG - 1), g);
      }
    }
  }
  for (int g = 0; g < ss_d.NS; ++g) {
    for (int hm = 0; hm < ss_c.hc2DS; ++hm) {
      for (int hu = 0; hu < ss_h.hTS; ++hu) {
        i_hc.age15_art_pop(hu, hm, g) += c_hc.hc2_art_pop(hu, hm, (ss_c.hc2AG - 1), g);
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_hiv_and_art_stratified_ageing(int t,
                                       const Parameters<ModelVariant, real_type> &pars,
                                       const State<ModelVariant, real_type> &state_curr,
                                       State<ModelVariant, real_type> &state_next,
                                       IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>();
  constexpr auto ss_h = ss.hiv;
  constexpr auto ss_d = ss.dp;
  constexpr auto ss_c = ss.children;
  const auto& p_op = pars.options;
  const auto& c_ha = state_curr.hiv;
  auto& n_ha = state_next.hiv;
  auto& i_ha = intermediate.hiv;
  auto& i_hc = intermediate.children;

  // age coarse stratified HIV population
  for (int g = 0; g < ss_d.NS; ++g) {
    int a = p_op.p_idx_hiv_first_adult;
    // Note: loop stops at hAG-1; no one ages out of the open-ended
    // age group
    for (int ha = 0; ha < (ss_h.hAG - 1); ++ha) {
      for (int i = 0; i < ss_h.hAG_span[ha]; ++i, ++a) {
        i_ha.hiv_age_up_prob(ha, g) += c_ha.p_hiv_pop(a, g);
      }

      if (i_ha.hiv_age_up_prob(ha, g) > 0) {
        i_ha.hiv_age_up_prob(ha, g) = c_ha.p_hiv_pop(a - 1, g) / i_ha.hiv_age_up_prob(ha, g);
      } else {
        i_ha.hiv_age_up_prob(ha, g) = 0.0;
      }
    }
  }


  for (int g = 0; g < ss_d.NS; ++g) {
    for (int ha = 1; ha < ss_h.hAG; ++ha) {
      for (int hm = 0; hm < ss_h.hDS; ++hm) {
        n_ha.h_hiv_adult(hm, ha, g) = ((1.0 - i_ha.hiv_age_up_prob(ha, g)) * c_ha.h_hiv_adult(hm, ha, g)) +
                                      (i_ha.hiv_age_up_prob(ha - 1, g) * c_ha.h_hiv_adult(hm, ha - 1, g));
        if (t > p_op.ts_art_start)
          for (int hu = 0; hu < ss_h.hTS; ++hu) {
            n_ha.h_art_adult(hu, hm, ha, g) = ((1.0 - i_ha.hiv_age_up_prob(ha, g)) * c_ha.h_art_adult(hu, hm, ha, g)) +
                                              (i_ha.hiv_age_up_prob(ha - 1, g) * c_ha.h_art_adult(hu, hm, ha - 1, g));
          }
      }
    }
  }

  if constexpr (ModelVariant::run_child_model) {
    const auto& p_hc = pars.children.children;

    for (int g = 0; g < ss_d.NS; ++g) {
      for (int hm = 0; hm < ss_h.hDS; ++hm) {
        for (int hm_adol = 0; hm_adol < ss_c.hc2DS; ++hm_adol){
          n_ha.h_hiv_adult(hm, 0, g) += i_hc.age15_hiv_pop(hm_adol, g) * p_hc.adult_cd4_dist(hm, hm_adol);
          if ((t > p_op.ts_art_start)) {
            for (int hu = 0; hu < ss_h.hTS; ++hu) {
              n_ha.h_art_adult(hu,hm, 0, g) += i_hc.age15_art_pop(hu, hm_adol, g) * p_hc.adult_cd4_dist(hm, hm_adol);
            }
          }
        }
      }
    }
  } else {
    for (int g = 0; g < ss_d.NS; ++g) {
      for (int hm = 0; hm < ss_h.hDS; ++hm) {
        n_ha.h_hiv_adult(hm, 0, g) = (1.0 - i_ha.hiv_age_up_prob(0, g)) * c_ha.h_hiv_adult(hm, 0, g);
        if (t > p_op.ts_art_start) {
          for (int hu = 0; hu < ss_h.hTS; ++hu) {
            n_ha.h_art_adult(hu, hm, 0, g) = (1.0 - i_ha.hiv_age_up_prob(0, g)) * c_ha.h_art_adult(hu, hm, 0, g);
          }
        }
      }
    }
  }
}


template<typename ModelVariant, typename real_type>
void run_hiv_and_art_stratified_deaths_and_migration(int t,
                                                     const Parameters<ModelVariant, real_type> &pars,
                                                     const State<ModelVariant, real_type> &state_curr,
                                                     State<ModelVariant, real_type> &state_next,
                                                     IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& i_ha = intermediate.hiv;
  auto& i_dp = intermediate.dp;

  for (int g = 0; g < ss_d.NS; ++g) {
    int a = p_op.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss_h.hAG; ++ha) {
      for (int i = 0; i < ss_h.hAG_span[ha]; ++i, ++a) {
        i_ha.p_hiv_pop_coarse_ages(ha, g) += n_ha.p_hiv_pop(a, g);
      }
    }
  }

  // remove non-HIV deaths and net migration from hiv stratified population
  for (int g = 0; g < ss_d.NS; ++g) {
    for (int a = 1; a < ss_d.pAG; ++a) {
      n_ha.p_hiv_pop(a, g) -= n_ha.p_hiv_pop_natural_deaths(a, g);
      if (p_op.proj_period_int == internal::PROJPERIOD_MIDYEAR) {
	      i_ha.hiv_net_migration(a, g) = n_ha.p_hiv_pop(a, g) * i_dp.migration_rate(a, g);
	      n_ha.p_hiv_pop(a, g) += i_ha.hiv_net_migration(a, g);
      }
    }
  }

  // remove non-HIV deaths and net migration from adult stratified population
  for (int g = 0; g < ss_d.NS; ++g) {
    int a = p_op.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss_h.hAG; ++ha) {
      real_type deaths_migrate = 0.0;
      for (int i = 0; i < ss_h.hAG_span[ha]; ++i, ++a) {
        deaths_migrate -= n_ha.p_hiv_pop_natural_deaths(a, g);
        if (p_op.proj_period_int == internal::PROJPERIOD_MIDYEAR) {
          deaths_migrate += i_ha.hiv_net_migration(a, g);
        }
      }

      real_type deaths_migrate_rate = 0.0;
      if (i_ha.p_hiv_pop_coarse_ages(ha, g) > 0) {
        deaths_migrate_rate = deaths_migrate / i_ha.p_hiv_pop_coarse_ages(ha, g);
      }

      for (int hm = 0; hm < ss_h.hDS; ++hm) {
        n_ha.h_hiv_adult(hm, ha, g) *= 1.0 + deaths_migrate_rate;
        if (t > p_op.ts_art_start) {
          for (int hu = 0; hu < ss_h.hTS; ++hu) {
            n_ha.h_art_adult(hu, hm, ha, g) *= 1.0 + deaths_migrate_rate;
          }
        }
      }
    }
  }
}

}

template<typename ModelVariant, typename real_type>
void run_hiv_pop_demographic_projection(int t,
                                        const Parameters<ModelVariant, real_type> &pars,
                                        const State<ModelVariant, real_type> &state_curr,
                                        State<ModelVariant, real_type> &state_next,
                                        internal::IntermediateData<ModelVariant, real_type> &intermediate) {

  internal::run_hiv_ageing_and_mortality<ModelVariant>(t, pars, state_curr, state_next, intermediate);

  if constexpr (ModelVariant::run_child_model) {
    internal::run_age_15_entrants<ModelVariant>(t, pars, state_curr, state_next, intermediate);
  }

  internal::run_hiv_and_art_stratified_ageing<ModelVariant>(t, pars, state_curr, state_next, intermediate);

  internal::run_hiv_and_art_stratified_deaths_and_migration<ModelVariant>(t, pars, state_curr, state_next,
                                                                          intermediate);
}


template<typename ModelVariant, typename real_type>
void run_hiv_pop_end_year_migration(int t,
                                    const Parameters<ModelVariant, real_type> &pars,
                                    const State<ModelVariant, real_type> &state_curr,
                                    State<ModelVariant, real_type> &state_next,
                                    internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& i_ha = intermediate.hiv;
  auto& i_dp = intermediate.dp;

  // remove net migration from hiv stratified population
  for (int g = 0; g < ss_d.NS; ++g) {
    for (int a = 0; a < ss_d.pAG; ++a) {
      i_ha.hiv_net_migration(a, g) = n_ha.p_hiv_pop(a, g) * i_dp.migration_rate(a, g);
      n_ha.p_hiv_pop(a, g) += i_ha.hiv_net_migration(a, g);
    }
  }

  // remove net migration from adult stratified population
  for (int g = 0; g < ss_d.NS; ++g) {
    int a = p_op.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss_h.hAG; ++ha) {
      real_type migration_num_ha = 0.0;
      real_type hivpop_ha_postmig = 0.0;
      for (int i = 0; i < ss_h.hAG_span[ha]; ++i, ++a) {
        hivpop_ha_postmig += n_ha.p_hiv_pop(a, g);
        migration_num_ha += i_ha.hiv_net_migration(a, g);
      }

      real_type migration_rate = 0.0;
      if (hivpop_ha_postmig > 0.0) {
        migration_rate = migration_num_ha / (hivpop_ha_postmig - migration_num_ha);
      }

      for (int hm = 0; hm < ss_h.hDS; ++hm) {
        n_ha.h_hiv_adult(hm, ha, g) *= 1.0 + migration_rate;
        if (t >= p_op.ts_art_start) {
          for (int hu = 0; hu < ss_h.hTS; ++hu) {
            n_ha.h_art_adult(hu, hm, ha, g) *= 1.0 + migration_rate;
          }
        }
      }
    }
  }
}

}
