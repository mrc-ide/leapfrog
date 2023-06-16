#pragma once

#include <iostream>
#include "types.hpp"

namespace leapfrog {

template<typename real_type>
void run_hiv_pop_demographic_projection(int time_step,
                                        const Parameters<real_type> &pars,
                                        const State<real_type> &state_curr,
                                        State<real_type> &state_next,
                                        internal::IntermediateData<real_type> &intermediate) {
  run_hiv_ageing_and_mortality(time_step, pars, state_curr, state_next,
                               intermediate);
  run_hiv_and_art_stratified_ageing(time_step, pars, state_curr, state_next,
                                    intermediate);
  run_hiv_and_art_stratified_deaths_and_migration(time_step, pars, state_curr,
                                                  state_next, intermediate);
}

namespace internal {

template<typename real_type>
void run_hiv_ageing_and_mortality(int time_step,
                                  const Parameters<real_type> &pars,
                                  const State<real_type> &state_curr,
                                  State<real_type> &state_next,
                                  IntermediateData<real_type> &intermediate) {
  // Non-hiv deaths
  for (int g = 0; g < pars.num_genders; ++g) {
    for (int a = 1; a < pars.age_groups_pop; ++a) {
      state_next.hiv_natural_deaths(a, g) =
          state_curr.hiv_population(a - 1, g) *
          (1.0 - pars.survival(a, g, time_step));
      state_next.hiv_population(a, g) = state_curr.hiv_population(a - 1, g);
    }

    // open age group
    state_next.hiv_natural_deaths(pars.age_groups_pop - 1, g) +=
        state_curr.hiv_population(pars.age_groups_pop - 1, g) *
        (1.0 - pars.survival(pars.age_groups_pop, g, time_step));
    state_next.hiv_population(pars.age_groups_pop - 1, g) +=
        state_curr.hiv_population(pars.age_groups_pop - 1, g);
  }
}

template<typename real_type>
void run_hiv_and_art_stratified_ageing(int time_step,
                                       const Parameters<real_type> &pars,
                                       const State<real_type> &state_curr,
                                       State<real_type> &state_next,
                                       IntermediateData<real_type> &intermediate) {
  // age coarse stratified HIV population
  for (int g = 0; g < pars.num_genders; ++g) {
    int a = pars.hiv_adult_first_age_group;
    // Note: loop stops at age_groups_hiv-1; no one ages out of the open-ended
    // age group
    for (int ha = 0; ha < (pars.age_groups_hiv - 1); ++ha) {
      for (int i = 0; i < pars.age_groups_hiv_span(ha); ++i, ++a) {
        intermediate.hiv_age_up_prob(ha, g) += state_curr.hiv_population(a, g);
      }
      if (intermediate.hiv_age_up_prob(ha, g) > 0) {
        intermediate.hiv_age_up_prob(ha, g) =
            state_curr.hiv_population(a - 1, g) / intermediate.hiv_age_up_prob(ha, g);
      } else {
        intermediate.hiv_age_up_prob(ha, g) = 0.0;
      }
    }
  }

  for (int g = 0; g < pars.num_genders; ++g) {
    for (int ha = 1; ha < pars.age_groups_hiv; ++ha) {
      for (int hm = 0; hm < pars.disease_stages; ++hm) {
        state_next.hiv_strat_adult(hm, ha, g) =
            ((1.0 - intermediate.hiv_age_up_prob(ha, g)) * state_curr.hiv_strat_adult(hm, ha, g)) +
            (intermediate.hiv_age_up_prob(ha - 1, g) * state_curr.hiv_strat_adult(hm, ha - 1, g));
        if (time_step > pars.time_art_start)
          for (int hu = 0; hu < pars.treatment_stages; ++hu) {
            state_next.art_strat_adult(hu, hm, ha, g) =
                ((1.0 - intermediate.hiv_age_up_prob(ha, g)) *
                 state_curr.art_strat_adult(hu, hm, ha, g)) +
                (intermediate.hiv_age_up_prob(ha - 1, g) *
                 state_curr.art_strat_adult(hu, hm, ha - 1, g));
          }
      }
    }
  }

  // TODO: add HIV+ 15 year old entrants see https://github.com/mrc-ide/leapfrog/issues/8
  for (int g = 0; g < pars.num_genders; ++g) {
    for (int hm = 0; hm < pars.disease_stages; ++hm) {
      state_next.hiv_strat_adult(hm, 0, g) =
          (1.0 - intermediate.hiv_age_up_prob(0, g)) * state_curr.hiv_strat_adult(hm, 0, g);
      // ADD HIV+ entrants here
      if (time_step > pars.time_art_start) {
        for (int hu = 0; hu < pars.treatment_stages; ++hu) {
          state_next.art_strat_adult(hu, hm, 0, g) =
              (1.0 - intermediate.hiv_age_up_prob(0, g)) *
              state_curr.art_strat_adult(hu, hm, 0, g);
          // ADD HIV+ entrants here
          //       artpop_t(hu, hm, 0, g, t) += paedsurv_g *
          //       paedsurv_artcd4dist(hu, hm, g, t) * entrantartcov(g, t);
        }
      }
    }
  }
}

template<typename real_type>
void run_hiv_and_art_stratified_deaths_and_migration(
    int time_step,
    const Parameters<real_type> &pars,
    const State<real_type> &state_curr,
    State<real_type> &state_next,
    IntermediateData<real_type> &intermediate) {
  for (int g = 0; g < pars.num_genders; ++g) {
    int a = pars.hiv_adult_first_age_group;
    for (int ha = 0; ha < pars.age_groups_hiv; ++ha) {
      for (int i = 0; i < pars.age_groups_hiv_span(ha); ++i, ++a) {
        intermediate.hiv_population_coarse_ages(ha, g) += state_next.hiv_population(a, g);
      }
    }
  }

  // remove non-HIV deaths and net migration from hiv stratified population
  for (int g = 0; g < pars.num_genders; ++g) {
    for (int a = 1; a < pars.age_groups_pop; ++a) {
      state_next.hiv_population(a, g) -= state_next.hiv_natural_deaths(a, g);
      intermediate.hiv_net_migration(a, g) =
          state_next.hiv_population(a, g) * intermediate.migration_rate(a, g);
      state_next.hiv_population(a, g) += intermediate.hiv_net_migration(a, g);
    }
  }

  // remove non-HIV deaths and net migration from adult stratified population
  for (int g = 0; g < pars.num_genders; ++g) {
    int a = pars.hiv_adult_first_age_group;
    for (int ha = 0; ha < pars.age_groups_hiv; ++ha) {
      real_type deaths_migrate = 0;
      for (int i = 0; i < pars.age_groups_hiv_span(ha); ++i, ++a) {
        deaths_migrate += (intermediate.hiv_net_migration(a, g) - state_next.hiv_natural_deaths(a, g));
      }

      real_type deaths_migrate_rate = 0.0;
      if (intermediate.hiv_population_coarse_ages(ha, g) > 0) {
        deaths_migrate_rate = deaths_migrate / intermediate.hiv_population_coarse_ages(ha, g);
      }

      for (int hm = 0; hm < pars.disease_stages; ++hm) {
        state_next.hiv_strat_adult(hm, ha, g) *= 1.0 + deaths_migrate_rate;
        if (time_step > pars.time_art_start) {
          for (int hu = 0; hu < pars.treatment_stages; ++hu) {
            state_next.art_strat_adult(hu, hm, ha, g) *=
                1.0 + deaths_migrate_rate;
          }
        }
      }
    }
  }
}

}
}
