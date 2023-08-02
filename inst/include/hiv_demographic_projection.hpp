#pragma once

#include <iostream>
#include "types.hpp"

namespace leapfrog {

namespace internal {

template<HivAgeStratification S, typename real_type>
void run_hiv_ageing_and_mortality(int time_step,
                                  const Parameters<real_type> &pars,
                                  const State<S, real_type> &state_curr,
                                  State<S, real_type> &state_next,
                                  IntermediateData<S, real_type> &intermediate) {
  const auto demog = pars.demography;
  constexpr auto ss = StateSpace<S>();
  // Non-hiv deaths
  for (int g = 0; g < ss.NS; ++g) {
    for (int a = 1; a < ss.pAG; ++a) {
      state_next.p_hiv_pop_natural_deaths(a, g) =
          state_curr.p_hiv_pop(a - 1, g) *
          (1.0 - demog.survival_probability(a, g, time_step));
      state_next.p_hiv_pop(a, g) = state_curr.p_hiv_pop(a - 1, g);
    }

    // open age group
    state_next.p_hiv_pop_natural_deaths(ss.pAG - 1, g) +=
        state_curr.p_hiv_pop(ss.pAG - 1, g) *
        (1.0 - demog.survival_probability(ss.pAG, g, time_step));
    state_next.p_hiv_pop(ss.pAG - 1, g) +=
        state_curr.p_hiv_pop(ss.pAG - 1, g);
  }
}

template<HivAgeStratification S, typename real_type>
void run_age_15_entrants(int time_step,
                         const Parameters<real_type> &pars,
                         const State<S, real_type> &state_curr,
                         State<S, real_type> &state_next,
                         IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();

  //TO DO: add ART entrants here
  for (int g = 0; g < ss.NS; ++g) {
    for (int hm = 0; hm < ss.hDS; ++hm) {
      for (int htm = 0; htm < ss.hcTT; ++htm) {
        intermediate.age15_hiv_pop(hm, g) += state_curr.hc2_hiv_pop(hm, htm, ss.hc2AG, g);
      }
    }
  }
}

template<HivAgeStratification S, typename real_type>
void run_hiv_and_art_stratified_ageing(int time_step,
                                       const Parameters<real_type> &pars,
                                       const State<S, real_type> &state_curr,
                                       State<S, real_type> &state_next,
                                       IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  // age coarse stratified HIV population
  for (int g = 0; g < ss.NS; ++g) {
    int a = pars.options.p_idx_hiv_first_adult;
    // Note: loop stops at hAG-1; no one ages out of the open-ended
    // age group
    for (int ha = 0; ha < (ss.hAG - 1); ++ha) {
      for (int i = 0; i < ss.hAG_span[ha]; ++i, ++a) {
        intermediate.hiv_age_up_prob(ha, g) += state_curr.p_hiv_pop(a, g);
      }
      if (intermediate.hiv_age_up_prob(ha, g) > 0) {
        intermediate.hiv_age_up_prob(ha, g) =
            state_curr.p_hiv_pop(a - 1, g) / intermediate.hiv_age_up_prob(ha, g);
      } else {
        intermediate.hiv_age_up_prob(ha, g) = 0.0;
      }
    }
  }

  for (int g = 0; g < ss.NS; ++g) {
    for (int ha = 1; ha < ss.hAG; ++ha) {
      for (int hm = 0; hm < ss.hDS; ++hm) {
        state_next.h_hiv_adult(hm, ha, g) =
            ((1.0 - intermediate.hiv_age_up_prob(ha, g)) * state_curr.h_hiv_adult(hm, ha, g)) +
            (intermediate.hiv_age_up_prob(ha - 1, g) * state_curr.h_hiv_adult(hm, ha - 1, g));
        if (time_step > pars.options.ts_art_start)
          for (int hu = 0; hu < ss.hTS; ++hu) {
            state_next.h_art_adult(hu, hm, ha, g) =
                ((1.0 - intermediate.hiv_age_up_prob(ha, g)) *
                 state_curr.h_art_adult(hu, hm, ha, g)) +
                (intermediate.hiv_age_up_prob(ha - 1, g) *
                 state_curr.h_art_adult(hu, hm, ha - 1, g));
          }
      }
    }
  }

  // TODO: add HIV+ 15 year old entrants see https://github.com/mrc-ide/leapfrog/issues/8
  for (int g = 0; g < ss.NS; ++g) {
    for (int hm = 0; hm < ss.hDS; ++hm) {
      state_next.h_hiv_adult(hm, 0, g) =
          (1.0 - intermediate.hiv_age_up_prob(0, g)) * state_curr.h_hiv_adult(hm, 0, g);
      // ADD HIV+ entrants here
      if (time_step > pars.options.ts_art_start) {
        for (int hu = 0; hu < ss.hTS; ++hu) {
          state_next.h_art_adult(hu, hm, 0, g) =
              (1.0 - intermediate.hiv_age_up_prob(0, g)) *
              state_curr.h_art_adult(hu, hm, 0, g);
          // ADD HIV+ entrants here
          //       artpop_t(hu, hm, 0, g, t) += paedsurv_g *
          //       paedsurv_artcd4dist(hu, hm, g, t) * entrantartcov(g, t);
        }
      }
    }
  }
}



template<HivAgeStratification S, typename real_type>
void run_hiv_and_art_stratified_deaths_and_migration(
    int time_step,
    const Parameters<real_type> &pars,
    const State<S, real_type> &state_curr,
    State<S, real_type> &state_next,
    IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  for (int g = 0; g < ss.NS; ++g) {
    int a = pars.options.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss.hAG; ++ha) {
      for (int i = 0; i < ss.hAG_span[ha]; ++i, ++a) {
        intermediate.p_hiv_pop_coarse_ages(ha, g) += state_next.p_hiv_pop(a, g);
      }
    }
  }

  // remove non-HIV deaths and net migration from hiv stratified population
  for (int g = 0; g < ss.NS; ++g) {
    for (int a = 1; a < ss.pAG; ++a) {
      state_next.p_hiv_pop(a, g) -= state_next.p_hiv_pop_natural_deaths(a, g);
      intermediate.hiv_net_migration(a, g) =
          state_next.p_hiv_pop(a, g) * intermediate.migration_rate(a, g);
      state_next.p_hiv_pop(a, g) += intermediate.hiv_net_migration(a, g);
    }
  }

  // remove non-HIV deaths and net migration from adult stratified population
  for (int g = 0; g < ss.NS; ++g) {
    int a = pars.options.p_idx_hiv_first_adult;
    for (int ha = 0; ha < ss.hAG; ++ha) {
      real_type deaths_migrate = 0;
      for (int i = 0; i < ss.hAG_span[ha]; ++i, ++a) {
        deaths_migrate += (intermediate.hiv_net_migration(a, g) - state_next.p_hiv_pop_natural_deaths(a, g));
      }

      real_type deaths_migrate_rate = 0.0;
      if (intermediate.p_hiv_pop_coarse_ages(ha, g) > 0) {
        deaths_migrate_rate = deaths_migrate / intermediate.p_hiv_pop_coarse_ages(ha, g);
      }

      for (int hm = 0; hm < ss.hDS; ++hm) {
        state_next.h_hiv_adult(hm, ha, g) *= 1.0 + deaths_migrate_rate;
        if (time_step > pars.options.ts_art_start) {
          for (int hu = 0; hu < ss.hTS; ++hu) {
            state_next.h_art_adult(hu, hm, ha, g) *=
                1.0 + deaths_migrate_rate;
          }
        }
      }
    }
  }
}

}

template<HivAgeStratification S, typename real_type>
void run_hiv_pop_demographic_projection(int time_step,
                                        const Parameters<real_type> &pars,
                                        const State<S, real_type> &state_curr,
                                        State<S, real_type> &state_next,
                                        internal::IntermediateData<S, real_type> &intermediate) {
  internal::run_hiv_ageing_and_mortality<S>(time_step, pars, state_curr, state_next,
                                            intermediate);
  if (pars.options.run_child_model) {
    internal::run_age_15_entrants(time_step, pars, state_curr, state_next, intermediate);
  }

  internal::run_hiv_and_art_stratified_ageing<S>(time_step, pars, state_curr, state_next,
                                                 intermediate);
  internal::run_hiv_and_art_stratified_deaths_and_migration<S>(time_step, pars, state_curr,
                                                               state_next, intermediate);
}

}
