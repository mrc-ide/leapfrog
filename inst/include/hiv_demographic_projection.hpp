#pragma once

#include <iostream>
#include "types.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void run_hiv_ageing_and_mortality(int time_step,
                                  const Parameters<ModelVariant, real_type> &pars,
                                  const State<ModelVariant, real_type> &state_curr,
                                  State<ModelVariant, real_type> &state_next,
                                  IntermediateData<ModelVariant, real_type> &intermediate) {
  const auto demog = pars.base.demography;
  constexpr auto ss = StateSpace<ModelVariant>().base;


  // Non-hiv deaths
  for (int g = 0; g < ss.NS; ++g) {
    for (int a = 1; a < ss.pAG; ++a) {
      state_next.base.p_hiv_pop_natural_deaths(a, g) =
          state_curr.base.p_hiv_pop(a - 1, g) *
          (1.0 - demog.survival_probability(a, g, time_step));
      state_next.base.p_hiv_pop(a, g) = state_curr.base.p_hiv_pop(a - 1, g);
    }

    // open age group
    state_next.base.p_hiv_pop_natural_deaths(ss.pAG - 1, g) +=
        state_curr.base.p_hiv_pop(ss.pAG - 1, g) *
        (1.0 - demog.survival_probability(ss.pAG, g, time_step));
    state_next.base.p_hiv_pop(ss.pAG - 1, g) +=
        state_curr.base.p_hiv_pop(ss.pAG - 1, g);
  }


}

template<typename ModelVariant, typename real_type>
void run_age_15_entrants(int time_step,
                         const Parameters<ModelVariant, real_type> &pars,
                         const State<ModelVariant, real_type> &state_curr,
                         State<ModelVariant, real_type> &state_next,
                         IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;

  for (int g = 0; g < ss.NS; ++g) {
    for (int hm = 0; hm < ss.hDS; ++hm) {
      for (int htm = 0; htm < hc_ss.hcTT; ++htm) {
        intermediate.children.age15_hiv_pop(hm, g) += state_curr.children.hc2_hiv_pop(hm, htm, (hc_ss.hc2AG-1), g);
      }
    }
  }
  for (int g = 0; g < ss.NS; ++g) {
    for (int hm = 0; hm < hc_ss.hc2DS; ++hm) {
      for (int hu = 0; hu < ss.hTS; ++hu) {
        intermediate.children.age15_art_pop(hu, hm, g) += state_curr.children.hc2_art_pop(hu, hm, (hc_ss.hc2AG-1), g);
      }
    }
  }

}

template<typename ModelVariant, typename real_type>
void run_hiv_and_art_stratified_ageing(int time_step,
                                       const Parameters<ModelVariant, real_type> &pars,
                                       const State<ModelVariant, real_type> &state_curr,
                                       State<ModelVariant, real_type> &state_next,
                                       IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children;
  // age coarse stratified HIV population
  for (int g = 0; g < ss.NS; ++g) {
    int a = pars.base.options.p_idx_hiv_first_adult;
    // Note: loop stops at hAG-1; no one ages out of the open-ended
    // age group
    for (int ha = 0; ha < (ss.hAG - 1); ++ha) {
      for (int i = 0; i < ss.hAG_span[ha]; ++i, ++a) {
        intermediate.base.hiv_age_up_prob(ha, g) += state_curr.base.p_hiv_pop(a, g);
      }
      if (intermediate.base.hiv_age_up_prob(ha, g) > 0) {
        intermediate.base.hiv_age_up_prob(ha, g) =
            state_curr.base.p_hiv_pop(a - 1, g) / intermediate.base.hiv_age_up_prob(ha, g);
      } else {
        intermediate.base.hiv_age_up_prob(ha, g) = 0.0;
      }
    }
  }


  for (int g = 0; g < ss.NS; ++g) {
    for (int ha = 1; ha < ss.hAG; ++ha) {
      for (int hm = 0; hm < ss.hDS; ++hm) {
        state_next.base.h_hiv_adult(hm, ha, g) =
            ((1.0 - intermediate.base.hiv_age_up_prob(ha, g)) * state_curr.base.h_hiv_adult(hm, ha, g)) +
            (intermediate.base.hiv_age_up_prob(ha - 1, g) * state_curr.base.h_hiv_adult(hm, ha - 1, g));
        if (time_step > pars.base.options.ts_art_start)
          for (int hu = 0; hu < ss.hTS; ++hu) {
            state_next.base.h_art_adult(hu, hm, ha, g) =
                ((1.0 - intermediate.base.hiv_age_up_prob(ha, g)) *
                 state_curr.base.h_art_adult(hu, hm, ha, g)) +
                (intermediate.base.hiv_age_up_prob(ha - 1, g) *
                 state_curr.base.h_art_adult(hu, hm, ha - 1, g));
          }
      }
    }
  }


  // TODO: add HIV+ 15 year old entrants see https://github.com/mrc-ide/leapfrog/issues/8
  if constexpr (ModelVariant::run_child_model) {
    for (int g = 0; g < ss.NS; ++g) {
      for (int hm = 0; hm < ss.hDS; ++hm) {
        for (int hm_adol = 0; hm_adol < hc_ss.hc2DS; ++hm_adol){
          state_next.base.h_hiv_adult(hm, 0, g) += intermediate.children.age15_hiv_pop(hm_adol, g) * cpars.children.adult_cd4_dist(hm, hm_adol);
          if ((time_step > pars.base.options.ts_art_start)) {
            for (int hu = 0; hu < ss.hTS; ++hu) {
              state_next.base.h_art_adult(hu,hm, 0, g) += intermediate.children.age15_art_pop(hu, hm_adol, g) * cpars.children.adult_cd4_dist(hm, hm_adol);
            }
          }
        }
      }
    }



  }else{
    for (int g = 0; g < ss.NS; ++g) {
      for (int hm = 0; hm < ss.hDS; ++hm) {
        state_next.base.h_hiv_adult(hm, 0, g) =
          (1.0 - intermediate.base.hiv_age_up_prob(0, g)) * state_curr.base.h_hiv_adult(hm, 0, g);
        if (time_step > pars.base.options.ts_art_start) {
          for (int hu = 0; hu < ss.hTS; ++hu) {
            state_next.base.h_art_adult(hu, hm, 0, g) = (1.0 - intermediate.base.hiv_age_up_prob(0, g)) *
              state_curr.base.h_art_adult(hu, hm, 0, g);
          }
        }
      }
    }
  }


  }


template<typename ModelVariant, typename real_type>
void run_hiv_and_art_stratified_deaths_and_migration(
    int time_step,
    const Parameters<ModelVariant, real_type> &pars,
    const State<ModelVariant, real_type> &state_curr,
    State<ModelVariant, real_type> &state_next,
    IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  for (int g = 0; g < ss.NS; ++g) {
    int a = pars.base.options.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss.hAG; ++ha) {
      for (int i = 0; i < ss.hAG_span[ha]; ++i, ++a) {
        intermediate.base.p_hiv_pop_coarse_ages(ha, g) += state_next.base.p_hiv_pop(a, g);
      }
    }
  }

  // remove non-HIV deaths and net migration from hiv stratified population
  for (int g = 0; g < ss.NS; ++g) {
    for (int a = 1; a < ss.pAG; ++a) {
      state_next.base.p_hiv_pop(a, g) -= state_next.base.p_hiv_pop_natural_deaths(a, g);
      intermediate.base.hiv_net_migration(a, g) =
          state_next.base.p_hiv_pop(a, g) * intermediate.base.migration_rate(a, g);
      state_next.base.p_hiv_pop(a, g) += intermediate.base.hiv_net_migration(a, g);
    }
  }

  // remove non-HIV deaths and net migration from adult stratified population
  for (int g = 0; g < ss.NS; ++g) {
    int a = pars.base.options.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss.hAG; ++ha) {
      real_type deaths_migrate = 0;
      for (int i = 0; i < ss.hAG_span[ha]; ++i, ++a) {
        deaths_migrate += (intermediate.base.hiv_net_migration(a, g) - state_next.base.p_hiv_pop_natural_deaths(a, g));
      }

      real_type deaths_migrate_rate = 0.0;
      if (intermediate.base.p_hiv_pop_coarse_ages(ha, g) > 0) {
        deaths_migrate_rate = deaths_migrate / intermediate.base.p_hiv_pop_coarse_ages(ha, g);
      }

      for (int hm = 0; hm < ss.hDS; ++hm) {
        state_next.base.h_hiv_adult(hm, ha, g) *= 1.0 + deaths_migrate_rate;
        if (time_step > pars.base.options.ts_art_start) {
          for (int hu = 0; hu < ss.hTS; ++hu) {
            state_next.base.h_art_adult(hu, hm, ha, g) *=
                1.0 + deaths_migrate_rate;
          }
        }
      }
    }
  }
}

}

template<typename ModelVariant, typename real_type>
void run_hiv_pop_demographic_projection(int time_step,
                                        const Parameters<ModelVariant, real_type> &pars,
                                        const State<ModelVariant, real_type> &state_curr,
                                        State<ModelVariant, real_type> &state_next,
                                        internal::IntermediateData<ModelVariant, real_type> &intermediate) {

  internal::run_hiv_ageing_and_mortality<ModelVariant>(time_step, pars, state_curr, state_next,
                                                       intermediate);

  if constexpr (ModelVariant::run_child_model) {
    internal::run_age_15_entrants<ModelVariant>(time_step, pars, state_curr, state_next, intermediate);
  }

  internal::run_hiv_and_art_stratified_ageing<ModelVariant>(time_step, pars, state_curr, state_next,
                                                            intermediate);

  internal::run_hiv_and_art_stratified_deaths_and_migration<ModelVariant>(time_step, pars, state_curr,
                                                                          state_next, intermediate);

}

}
