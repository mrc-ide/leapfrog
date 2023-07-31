#pragma once

#include "types.hpp"

namespace leapfrog {

template<HivAgeStratification S, typename real_type>
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
    Tensor5<real_type> hc1_hiv_pop;
    Tensor5<real_type> hc2_hiv_pop;
    Tensor5<real_type> hc1_art_pop;
    Tensor5<real_type> hc2_art_pop;
    Tensor5<real_type> hc1_art_aids_deaths;
    Tensor5<real_type> hc2_art_aids_deaths;
    Tensor5<real_type> hc1_noart_aids_deaths;
    Tensor5<real_type> hc2_noart_aids_deaths;

    OutputState(int no_output_years)
        : total_population(StateSpace<S>().age_groups_pop, StateSpace<S>().num_genders, no_output_years),
          natural_deaths(StateSpace<S>().age_groups_pop, StateSpace<S>().num_genders, no_output_years),
          hiv_population(StateSpace<S>().age_groups_pop, StateSpace<S>().num_genders, no_output_years),
          hiv_natural_deaths(StateSpace<S>().age_groups_pop, StateSpace<S>().num_genders, no_output_years),
          hiv_strat_adult(StateSpace<S>().disease_stages, StateSpace<S>().age_groups_hiv, StateSpace<S>().num_genders,
                          no_output_years),
          art_strat_adult(StateSpace<S>().treatment_stages,
                          StateSpace<S>().disease_stages,
                          StateSpace<S>().age_groups_hiv,
                          StateSpace<S>().num_genders,
                          no_output_years),
          births(no_output_years),
          aids_deaths_no_art(StateSpace<S>().disease_stages, StateSpace<S>().age_groups_hiv,
                             StateSpace<S>().num_genders, no_output_years),
          infections(StateSpace<S>().age_groups_pop, StateSpace<S>().num_genders, no_output_years),
          aids_deaths_art(StateSpace<S>().treatment_stages, StateSpace<S>().disease_stages,
                          StateSpace<S>().age_groups_hiv, StateSpace<S>().num_genders, no_output_years),
          art_initiation(StateSpace<S>().disease_stages, StateSpace<S>().age_groups_hiv,
                         StateSpace<S>().num_genders, no_output_years),
          hiv_deaths(StateSpace<S>().age_groups_pop, StateSpace<S>().num_genders, no_output_years),
          hc1_hiv_pop(StateSpace<S>().hc1_disease_stages, StateSpace<S>().hTM, StateSpace<S>().hc1_age_groups,
                     StateSpace<S>().num_genders, no_output_years),
          hc2_hiv_pop(StateSpace<S>().hc2_disease_stages, StateSpace<S>().hTM, StateSpace<S>().hc2_age_groups,
                     StateSpace<S>().num_genders, no_output_years),
          hc1_art_pop(StateSpace<S>().treatment_stages, StateSpace<S>().hc1_disease_stages, StateSpace<S>().hc1_age_groups,
                                 StateSpace<S>().num_genders, no_output_years),
          hc2_art_pop(StateSpace<S>().treatment_stages, StateSpace<S>().hc2_disease_stages, StateSpace<S>().hc2_age_groups,
                                             StateSpace<S>().num_genders, no_output_years),
          hc1_art_aids_deaths(StateSpace<S>().treatment_stages, StateSpace<S>().hc1_disease_stages, StateSpace<S>().hc1_age_groups,
                                             StateSpace<S>().num_genders, no_output_years),
          hc2_art_aids_deaths(StateSpace<S>().treatment_stages, StateSpace<S>().hc2_disease_stages, StateSpace<S>().hc2_age_groups,
                                   StateSpace<S>().num_genders, no_output_years),
          hc1_noart_aids_deaths(StateSpace<S>().hc1_disease_stages, StateSpace<S>().hTM, StateSpace<S>().hc1_age_groups,
                                   StateSpace<S>().num_genders, no_output_years),
          hc2_noart_aids_deaths(StateSpace<S>().hc2_disease_stages, StateSpace<S>().hTM, StateSpace<S>().hc2_age_groups,
                                                       StateSpace<S>().num_genders, no_output_years)
                                             {
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
      hc1_hiv_pop.setZero();
      hc2_hiv_pop.setZero();
      hc1_art_pop.setZero();
      hc2_art_pop.setZero();
      hc1_art_aids_deaths.setZero();
      hc2_art_aids_deaths.setZero();
      hc1_noart_aids_deaths.setZero();
      hc2_noart_aids_deaths.setZero();
    }
  };

  StateSaver(int time_steps,
             std::vector<int> save_steps) :
      save_steps(save_steps),
      full_state(save_steps.size()) {
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


  void save_state(const State<S, real_type> &state, int current_year) {
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
        full_state.hc1_hiv_pop.chip(i, full_state.hc1_hiv_pop.NumDimensions - 1) = state.hc1_hiv_pop;
        full_state.hc2_hiv_pop.chip(i, full_state.hc2_hiv_pop.NumDimensions - 1) = state.hc2_hiv_pop;
        full_state.hc1_art_pop.chip(i, full_state.hc1_art_pop.NumDimensions - 1) = state.hc1_art_pop;
        full_state.hc2_art_pop.chip(i, full_state.hc2_art_pop.NumDimensions - 1) = state.hc2_art_pop;
        full_state.hc1_art_aids_deaths.chip(i, full_state.hc1_art_aids_deaths.NumDimensions - 1) = state.hc1_art_aids_deaths;
        full_state.hc2_art_aids_deaths.chip(i, full_state.hc2_art_aids_deaths.NumDimensions - 1) = state.hc2_art_aids_deaths;
        full_state.hc1_noart_aids_deaths.chip(i, full_state.hc1_noart_aids_deaths.NumDimensions - 1) = state.hc1_noart_aids_deaths;
        full_state.hc2_noart_aids_deaths.chip(i, full_state.hc2_noart_aids_deaths.NumDimensions - 1) = state.hc2_noart_aids_deaths;
        return;
      }
    }
  }

  const OutputState &get_full_state() const {
    return full_state;
  }


private:
  std::vector<int> save_steps;
  OutputState full_state;
};

}
