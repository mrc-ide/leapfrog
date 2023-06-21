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
  };

  StateSaver(std::vector<int> save_steps,
             int age_groups_pop,
             int num_genders,
             int disease_stages,
             int age_groups_hiv,
             int treatment_stages) :
      save_steps(save_steps) {

    size_t no_output_years = save_steps.size();
    std::cout << "no output years " << no_output_years << "\n";
    std::cout << "initialising tot pop \n";
    full_state.total_population(age_groups_pop, num_genders, no_output_years);
    std::cout << "initialising nat deaths \n";
    full_state.natural_deaths(age_groups_pop, num_genders, no_output_years);
    std::cout << "initialising hiv pop \n";
    full_state.hiv_population(age_groups_pop, num_genders, no_output_years);
    std::cout << "initialising hiv nat deaths \n";
    full_state.hiv_natural_deaths(age_groups_pop, num_genders, no_output_years);
    std::cout << "initialising hiv strat adult \n";
    full_state.hiv_strat_adult(disease_stages, age_groups_hiv, num_genders, no_output_years);
    std::cout << "initialising art strat \n";
    full_state.art_strat_adult(treatment_stages,
                               disease_stages,
                               age_groups_hiv,
                               num_genders,
                               no_output_years);
    std::cout << "initialising births \n";
    full_state.births(no_output_years);
    std::cout << "initialising aids deaths no \n";
    full_state.aids_deaths_no_art(disease_stages, age_groups_hiv, num_genders, no_output_years);
    std::cout << "initialising ifnections \n";
    full_state.infections(age_groups_pop, num_genders, no_output_years);
    std::cout << "initialising aids deaths art \n";
    full_state.aids_deaths_art(treatment_stages, disease_stages, age_groups_hiv, num_genders, no_output_years);
    std::cout << "initialising art init \n";
    full_state.art_initiation(disease_stages, age_groups_hiv, num_genders, no_output_years);
    std::cout << "initialising hiv deaths \n";
    full_state.hiv_deaths(age_groups_pop, num_genders, no_output_years);
    std::cout << "setting 0 \n";

    full_state.total_population.setZero();
    full_state.natural_deaths.setZero();
    full_state.hiv_population.setZero();
    full_state.hiv_natural_deaths.setZero();
    full_state.hiv_strat_adult.setZero();
    full_state.art_strat_adult.setZero();
    full_state.births.setZero();
    full_state.aids_deaths_no_art.setZero();
    full_state.infections.setZero();
    full_state.aids_deaths_art.setZero();
    full_state.art_initiation.setZero();
    full_state.hiv_deaths.setZero();
  }


  void save_state(const State<real_type> &state, int current_year) {
    // Always report the first year and then every nth year after that
    for (int i = 0; i < save_steps.size(); ++i) {
      if (current_year == save_steps[i]) {
//        full_state.total_population.chip(i, 3) = state.total_population;
        std::cout << "reporting at time " << current_year << "\n";
        return;
      }
    }
//      full_state.natural_deaths[, , report_index] = state.natural_deaths;
//      full_state.hiv_population[, , report_index] = state.hiv_population;
//      full_state.hiv_natural_deaths[, , report_index] = state.hiv_natural_deaths;
//      full_state.hiv_strat_adult[, , , report_index] = state.hiv_strat_adult;
//      full_state.art_strat_adult[, , , , report_index] = state.art_strat_adult;
//      full_state.births[report_index] = state.births;
//      full_state.aids_deaths_no_art[, , , report_index] = state.aids_deaths_no_art;
//      full_state.infections[, , report_index] = state.infections;
//      full_state.aids_deaths_art[, , , , report_index] = state.aids_deaths_art;
//      full_state.art_initiation[, , , report_index] = state.art_initiation;
//      full_state.hiv_deaths[, , report_index] = state.hiv_deaths;
  }

  OutputState get_full_state() const {
    return full_state;
  }


private:
  std::vector<int> save_steps;
  OutputState full_state;
};

}
