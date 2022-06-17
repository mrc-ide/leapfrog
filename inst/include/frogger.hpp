#pragma once

#include "demographic_projection.hpp"

template <typename real_type>
State<real_type> run_model(int time_steps, const Parameters<real_type>& pars) {
  State<double> state(pars.age_groups_pop, pars.num_genders,
                      pars.disease_stages, pars.age_groups_hiv,
                      pars.treatment_stages);

  initialise_model_state(pars, state);
  auto state_next = state;
  WorkingData<real_type> working(pars.age_groups_pop, pars.num_genders);
  // Each time step is mid-point of the year
  for (int step = 1; step <= time_steps; ++step) {
    run_demographic_projection(step, pars, state, state_next, working);
    run_hiv_pop_demographic_projection(step, pars, state, state_next, working);
    std::swap(state, state_next);
    working.reset();
  }
  return state;
}

template <typename real_type>
void initialise_model_state(const Parameters<real_type>& pars,
                            State<real_type>& state) {
  for (int g = 0; g < pars.num_genders; g++) {
    for (int a = 0; a < pars.age_groups_pop; a++) {
      state.total_population(a, g) = pars.base_pop(a, g);
      state.natural_deaths(a, g) = 0;
    }
  }
  state.births = 0;
}

template <typename real_type>
void run_hiv_pop_demographic_projection(int time_step,
                                        const Parameters<real_type>& pars,
                                        const State<real_type>& state_curr,
                                        State<real_type>& state_next,
                                        WorkingData<real_type>& working) {
  const int num_genders = pars.num_genders;
  const int age_groups_pop = pars.age_groups_pop;
  const int age_groups_hiv = pars.age_groups_hiv;
  const int disease_stages = pars.disease_stages;
  const int treatment_stages = pars.treatment_stages;
  const int time_art_start = pars.time_art_start;

  // age population and calculate non-HIV deaths to HIV population
  for (int g = 0; g < num_genders; ++g) {
    for (int a = 1; a < age_groups_pop; ++a) {
      state_next.hiv_natural_deaths(a, g) =
          state_curr.hiv_population(a - 1, g) *
          (1.0 - pars.survival(a, g, time_step));
      state_next.hiv_population(a, g) = state_curr.hiv_population(a - 1, g);
    }

    // open age group
    state_next.hiv_natural_deaths(age_groups_pop - 1, g) +=
        state_curr.hiv_population(age_groups_pop - 1, g) *
        (1.0 - pars.survival(age_groups_pop, g, time_step));
    state_next.hiv_population(age_groups_pop - 1, g) +=
        state_curr.hiv_population(age_groups_pop - 1, g);
  }

  // age coarse stratified HIV population
  Tensor2<real_type> hiv_ag_prob(age_groups_hiv, num_genders);
  hiv_ag_prob.setZero();

  for (int g = 0; g < num_genders; ++g) {
    int a = pars.hiv_adult_first_age_group;
    // Note: loop stops at age_groups_hiv-1; no one ages out of the open-ended
    // age group
    for (int ha = 0; ha < (age_groups_hiv - 1); ++ha) {
      for (int i = 0; i < pars.age_groups_hiv_span(ha); ++i) {
        hiv_ag_prob(ha, g) += state_curr.hiv_population(a, g);
        a++;
      }
      hiv_ag_prob(ha, g) =
          (hiv_ag_prob(ha, g) > 0)
              ? state_curr.hiv_population(a - 1, g) / hiv_ag_prob(ha, g)
              : 0.0;
    }
  }

  for (int g = 0; g < num_genders; ++g) {
    for (int ha = 1; ha < age_groups_hiv; ++ha) {
      for (int hm = 0; hm < disease_stages; ++hm) {
        state_next.hiv_strat_adult(hm, ha, g) =
            (1.0 - hiv_ag_prob(ha, g)) *
            state_curr.hiv_strat_adult(hm, ha, g);  // age-out
        state_next.hiv_strat_adult(hm, ha, g) +=
            hiv_ag_prob(ha - 1, g) *
            state_curr.hiv_strat_adult(hm, ha - 1, g);  // age-in
        if (time_step > time_art_start)
          for (int hu = 0; hu < treatment_stages; hu++) {
            state_next.art_strat_adult(hu, hm, ha, g) =
                (1.0 - hiv_ag_prob(ha, g)) *
                state_curr.art_strat_adult(hu, hm, ha, g);
            state_next.art_strat_adult(hu, hm, ha, g) +=
                hiv_ag_prob(ha - 1, g) *
                state_curr.art_strat_adult(hu, hm, ha - 1, g);
          }
      }
    }
  }

  // TODO: add HIV+ 15 year old entrants
  for (int g = 0; g < num_genders; ++g) {
    for (int hm = 0; hm < disease_stages; ++hm) {
      state_next.hiv_strat_adult(hm, 0, g) =
          (1.0 - hiv_ag_prob(0, g)) * state_curr.hiv_strat_adult(hm, 0, g);
      // ADD HIV+ entrants here
      if (time_step > time_art_start) {
        for (int hu = 0; hu < treatment_stages; hu++) {
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

  Tensor2<real_type> hiv_population_ha(age_groups_hiv, num_genders);
  hiv_population_ha.setZero();
  for (int g = 0; g < num_genders; ++g) {
    int a = pars.hiv_adult_first_age_group;
    for (int ha = 0; ha < age_groups_hiv; ha++) {
      for (int i = 0; i < pars.age_groups_hiv_span(ha); i++) {
        hiv_population_ha(ha, g) += state_next.hiv_population(a, g);
        a++;
      }
    }
  }

  // remove non-HIV deaths and net migration from hiv_population
  Tensor2<real_type> net_migration_ag(age_groups_pop, num_genders);
  for (int g = 0; g < num_genders; g++) {
    for (int a = 1; a < age_groups_pop; a++) {
      state_next.hiv_population(a, g) -= state_next.hiv_natural_deaths(a, g);
      // net_migration_ag(a, g) =
      //     state_next.hiv_population(a, g) * migration_rate(a, g);
      state_next.hiv_population(a, g) += net_migration_ag(a, g);
    }
  }

  // remove non-HIV deaths and net migration from adult stratified population
  for (int g = 0; g < num_genders; g++) {
    int a = pars.hiv_adult_first_age_group;
    for (int ha = 0; ha < age_groups_hiv; ha++) {
      real_type deaths_migrate = 0;
      for (int i = 0; i < pars.age_groups_hiv_span(ha); i++) {
        deaths_migrate -= state_next.hiv_natural_deaths(a, g);
        deaths_migrate += net_migration_ag(a, g);
        a++;
      }

      real_type deaths_migrate_rate =
          hiv_population_ha(ha, g) > 0
              ? deaths_migrate / hiv_population_ha(ha, g)
              : 0.0;
      for (int hm = 0; hm < disease_stages; hm++) {
        state_next.hiv_strat_adult(hm, ha, g) *= 1.0 + deaths_migrate_rate;
        if (time_step > time_art_start) {
          for (int hu = 0; hu < treatment_stages; hu++) {
            state_next.art_strat_adult(hu, hm, ha, g) *=
                1.0 + deaths_migrate_rate;
          }
        }
      }
    }
  }
}
=======
>>>>>>> mrc-3273
