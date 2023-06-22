#pragma once

#include "types.hpp"

namespace leapfrog {

template<typename real_type>
class StateSaver {

public:

  struct OutputState {
    Tensor3<real_type> total_population;
    Tensor3<real_type> natural_deaths;
    Tensor3<real_type> hiv_population;
    Tensor3<real_type> hiv_natural_deaths;
    Tensor4<real_type> hiv_strat_adult;
    Tensor5<real_type> art_strat_adult;
    Tensor1<real_type> births;
    Tensor4<real_type> aids_deaths_no_art;
    Tensor3<real_type> infections;
    Tensor5<real_type> aids_deaths_art;
    Tensor4<real_type> art_initiation;
    Tensor3<real_type> hiv_deaths;

    OutputState(int age_groups_pop,
                int num_genders,
                int disease_stages,
                int age_groups_hiv,
                int treatment_stages,
                int no_output_years)
        : total_population(age_groups_pop, num_genders, no_output_years),
          natural_deaths(age_groups_pop, num_genders, no_output_years),
          hiv_population(age_groups_pop, num_genders, no_output_years),
          hiv_natural_deaths(age_groups_pop, num_genders, no_output_years),
          hiv_strat_adult(disease_stages, age_groups_hiv, num_genders, no_output_years),
          art_strat_adult(treatment_stages,
                          disease_stages,
                          age_groups_hiv,
                          num_genders,
                          no_output_years),
          births(no_output_years),
          aids_deaths_no_art(disease_stages, age_groups_hiv, num_genders, no_output_years),
          infections(age_groups_pop, num_genders, no_output_years),
          aids_deaths_art(treatment_stages, disease_stages, age_groups_hiv, num_genders, no_output_years),
          art_initiation(disease_stages, age_groups_hiv, num_genders, no_output_years),
          hiv_deaths(age_groups_pop, num_genders, no_output_years) {

      total_population.setZero();
      natural_deaths.setZero();
      hiv_population.setZero();
      hiv_natural_deaths.setZero();
      hiv_strat_adult.setZero();
      art_strat_adult.setZero();
      births.setZero();
      aids_deaths_no_art.setZero();
      infections.setZero();
      aids_deaths_art.setZero();
      art_initiation.setZero();
      hiv_deaths.setZero();
    }
  };

  StateSaver(int time_steps,
             std::vector<int> save_steps,
             int age_groups_pop,
             int num_genders,
             int disease_stages,
             int age_groups_hiv,
             int treatment_stages) :
      save_steps(save_steps),
      full_state(age_groups_pop, num_genders, disease_stages, age_groups_hiv, treatment_stages, save_steps.size()) {
    for (int step: save_steps) {
      if (step < 0) {
        std::stringstream ss;
        ss << "Output step must be at least 0, got '" << step << "'." << std::endl;
        throw std::runtime_error(ss.str());
      }
      if (step > time_steps) {
        std::stringstream ss;
        ss << "Output step can be at most number of time steps run which is '" << time_steps << "', got step '" << step
           << "'." << std::endl;
        throw std::runtime_error(ss.str());
      }
    }
  }


  void save_state(const State<real_type> &state, int current_year) {
    for (size_t i = 0; i < save_steps.size(); ++i) {
      if (current_year == save_steps[i]) {
        full_state.total_population.chip(i, full_state.total_population.NumDimensions - 1) = state.total_population;
        full_state.natural_deaths.chip(i, full_state.natural_deaths.NumDimensions - 1) = state.natural_deaths;
        full_state.hiv_population.chip(i, full_state.hiv_population.NumDimensions - 1) = state.hiv_population;
        full_state.hiv_natural_deaths.chip(i, full_state.hiv_natural_deaths.NumDimensions - 1) =
            state.hiv_natural_deaths;
        full_state.hiv_strat_adult.chip(i, full_state.hiv_strat_adult.NumDimensions - 1) = state.hiv_strat_adult;
        full_state.art_strat_adult.chip(i, full_state.art_strat_adult.NumDimensions - 1) = state.art_strat_adult;
        full_state.births(i) = state.births;
        full_state.aids_deaths_no_art.chip(i, full_state.aids_deaths_no_art.NumDimensions - 1) =
            state.aids_deaths_no_art;
        full_state.infections.chip(i, full_state.infections.NumDimensions - 1) = state.infections;
        full_state.aids_deaths_art.chip(i, full_state.aids_deaths_art.NumDimensions - 1) = state.aids_deaths_art;
        full_state.art_initiation.chip(i, full_state.art_initiation.NumDimensions - 1) = state.art_initiation;
        full_state.hiv_deaths.chip(i, full_state.hiv_deaths.NumDimensions - 1) = state.hiv_deaths;
        return;
      }
    }
  }

  OutputState get_full_state() const {
    return full_state;
  }


private:
  std::vector<int> save_steps;
  OutputState full_state;
};

}
