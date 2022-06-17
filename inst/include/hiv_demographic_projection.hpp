#pragma once

#include <iostream>
#include "types.hpp"

template <typename real_type>
void run_hiv_pop_demographic_projection(int time_step,
                                        const Parameters<real_type>& pars,
                                        const State<real_type>& state_curr,
                                        State<real_type>& state_next,
                                        WorkingData<real_type>& working) {
  run_hiv_ageing_and_mortality(time_step, pars, state_curr, state_next);
  run_hiv_and_art_stratified_ageing(time_step, pars, state_curr, state_next);
  run_hiv_deaths_and_migration(time_step, pars, state_next, working);
  run_hiv_and_art_stratified_deaths_and_migration(time_step, pars, state_curr,
                                                  state_next, working);
}

template <typename real_type>
void run_hiv_ageing_and_mortality(int time_step,
                                  const Parameters<real_type>& pars,
                                  const State<real_type>& state_curr,
                                  State<real_type>& state_next) {
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

template <typename real_type>
void run_hiv_and_art_stratified_ageing(int time_step,
                                       const Parameters<real_type>& pars,
                                       const State<real_type>& state_curr,
                                       State<real_type>& state_next) {
  // TODO: better name for hiv_ag_prob - what is this? Prob of chance of someone
  // getting HIV as a function of the size of the HIV population
  Tensor2<real_type> hiv_ag_prob(pars.age_groups_hiv, pars.num_genders);
  hiv_ag_prob.setZero();
  // age coarse stratified HIV population
  for (int g = 0; g < pars.num_genders; ++g) {
    int a = pars.hiv_adult_first_age_group;
    // Note: loop stops at age_groups_hiv-1; no one ages out of the open-ended
    // age group
    for (int ha = 0; ha < (pars.age_groups_hiv - 1); ++ha) {
      for (int i = 0; i < pars.age_groups_hiv_span(ha); ++i) {
        hiv_ag_prob(ha, g) += state_curr.hiv_population(a, g);
        a++;
      }
      hiv_ag_prob(ha, g) =
          (hiv_ag_prob(ha, g) > 0)
              ? state_curr.hiv_population(a - 1, g) / hiv_ag_prob(ha, g)
              : 0.0;
    }

    for (int ha = 1; ha < pars.age_groups_hiv; ++ha) {
      for (int hm = 0; hm < pars.disease_stages; ++hm) {
        state_next.hiv_strat_adult(hm, ha, g) =
            (1.0 - hiv_ag_prob(ha, g)) *
            state_curr.hiv_strat_adult(hm, ha, g);  // age-out
        state_next.hiv_strat_adult(hm, ha, g) +=
            hiv_ag_prob(ha - 1, g) *
            state_curr.hiv_strat_adult(hm, ha - 1, g);  // age-in
        if (time_step > pars.time_art_start)
          for (int hu = 0; hu < pars.treatment_stages; hu++) {
            state_next.art_strat_adult(hu, hm, ha, g) =
                (1.0 - hiv_ag_prob(ha, g)) *
                state_curr.art_strat_adult(hu, hm, ha, g);
            state_next.art_strat_adult(hu, hm, ha, g) +=
                hiv_ag_prob(ha - 1, g) *
                state_curr.art_strat_adult(hu, hm, ha - 1, g);
          }
      }
    }

    // TODO: add HIV+ 15 year old entrants
    for (int hm = 0; hm < pars.disease_stages; ++hm) {
      state_next.hiv_strat_adult(hm, 0, g) =
          (1.0 - hiv_ag_prob(0, g)) * state_curr.hiv_strat_adult(hm, 0, g);
      // ADD HIV+ entrants here
      if (time_step > pars.time_art_start) {
        for (int hu = 0; hu < pars.treatment_stages; hu++) {
          state_next.art_strat_adult(hu, hm, 0, g) =
              (1.0 - hiv_ag_prob(0, g)) *
              state_curr.art_strat_adult(hu, hm, 0, g);
          // ADD HIV+ entrants here
          //       artpop_t(hu, hm, 0, g, t) += paedsurv_g *
          //       paedsurv_artcd4dist(hu, hm, g, t) * entrantartcov(g, t);
        }
      }
    }
  }
}

template <typename real_type>
void run_hiv_deaths_and_migration(int time_step,
                                  const Parameters<real_type>& pars,
                                  State<real_type>& state_next,
                                  WorkingData<real_type>& working) {
  // remove non-HIV deaths and net migration from hiv_population
  for (int g = 0; g < pars.num_genders; g++) {
    for (int a = 1; a < pars.age_groups_pop; a++) {
      state_next.hiv_population(a, g) -= state_next.hiv_natural_deaths(a, g);
      working.hiv_net_migration(a, g) =
          state_next.hiv_population(a, g) * working.migration_rate(a, g);
      state_next.hiv_population(a, g) += working.hiv_net_migration(a, g);
    }
  }
}

template <typename real_type>
void run_hiv_and_art_stratified_deaths_and_migration(
    int time_step,
    const Parameters<real_type>& pars,
    const State<real_type>& state_curr,
    State<real_type>& state_next,
    WorkingData<real_type>& working) {
  // TODO: better name for hiv_population_ha?
  Tensor2<real_type> hiv_population_ha(pars.age_groups_hiv, pars.num_genders);
  hiv_population_ha.setZero();
  for (int g = 0; g < pars.num_genders; g++) {
    int a = pars.hiv_adult_first_age_group;
    for (int ha = 0; ha < pars.age_groups_hiv; ha++) {
      for (int i = 0; i < pars.age_groups_hiv_span(ha); i++) {
        hiv_population_ha(ha, g) += state_next.hiv_population(a, g);
        a++;
      }
    }
  }

  // remove non-HIV deaths and net migration from adult stratified population
  for (int g = 0; g < pars.num_genders; g++) {
    int a = pars.hiv_adult_first_age_group;
    for (int ha = 0; ha < pars.age_groups_hiv; ha++) {
      real_type deaths_migrate = 0;
      for (int i = 0; i < pars.age_groups_hiv_span(ha); i++) {
        deaths_migrate -= state_next.hiv_natural_deaths(a, g);
        deaths_migrate += working.hiv_net_migration(a, g);
        a++;
      }

      real_type deaths_migrate_rate =
          hiv_population_ha(ha, g) > 0
              ? deaths_migrate / hiv_population_ha(ha, g)
              : 0.0;
      for (int hm = 0; hm < pars.disease_stages; hm++) {
        state_next.hiv_strat_adult(hm, ha, g) *= 1.0 + deaths_migrate_rate;
        if (time_step > pars.time_art_start) {
          for (int hu = 0; hu < pars.treatment_stages; hu++) {
            state_next.art_strat_adult(hu, hm, ha, g) *=
                1.0 + deaths_migrate_rate;
          }
        }
      }
    }
  }
}
