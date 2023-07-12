#pragma once

#include <Rcpp.h>
#include <unsupported/Eigen/CXX11/Tensor>

#include "state_saver.hpp"

template<typename real_type, leapfrog::HivAgeStratification S>
Rcpp::List build_r_output(const typename leapfrog::StateSaver<real_type, S>::OutputState &state,
                          const std::vector<int> save_steps) {
  constexpr auto ss = leapfrog::StateSpace<S>();
  size_t output_years = save_steps.size();
  Rcpp::NumericVector r_total_population(ss.age_groups_pop * ss.num_genders * output_years);
  Rcpp::NumericVector r_births(output_years);
  Rcpp::NumericVector r_natural_deaths(ss.age_groups_pop * ss.num_genders * output_years);
  Rcpp::NumericVector r_hiv_population(ss.age_groups_pop * ss.num_genders * output_years);
  Rcpp::NumericVector r_hiv_natural_deaths(ss.age_groups_pop * ss.num_genders * output_years);
  Rcpp::NumericVector r_hiv_strat_adult(ss.disease_stages * ss.age_groups_hiv *
                                        ss.num_genders * output_years);
  Rcpp::NumericVector r_art_strat_adult(ss.treatment_stages * ss.disease_stages *
                                        ss.age_groups_hiv * ss.num_genders * output_years);
  Rcpp::NumericVector r_aids_deaths_no_art(ss.disease_stages * ss.age_groups_hiv * ss.num_genders * output_years);
  Rcpp::NumericVector r_infections(ss.age_groups_pop * ss.num_genders * output_years);
  Rcpp::NumericVector r_aids_deaths_art(ss.treatment_stages * ss.disease_stages * ss.age_groups_hiv *
                                        ss.num_genders * output_years);
  Rcpp::NumericVector r_art_initiation(ss.disease_stages * ss.age_groups_hiv * ss.num_genders * output_years);
  Rcpp::NumericVector r_hiv_deaths(ss.age_groups_pop * ss.num_genders * output_years);

  r_total_population.attr("dim") =
      Rcpp::NumericVector::create(ss.age_groups_pop, ss.num_genders, output_years);
  r_natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(ss.age_groups_pop, ss.num_genders, output_years);
  r_hiv_population.attr("dim") =
      Rcpp::NumericVector::create(ss.age_groups_pop, ss.num_genders, output_years);
  r_hiv_natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(ss.age_groups_pop, ss.num_genders, output_years);
  r_hiv_strat_adult.attr("dim") =
      Rcpp::NumericVector::create(ss.disease_stages, ss.age_groups_hiv, ss.num_genders, output_years);
  r_art_strat_adult.attr("dim") = Rcpp::NumericVector::create(
      ss.treatment_stages, ss.disease_stages, ss.age_groups_hiv, ss.num_genders, output_years);
  r_births.attr("dim") = Rcpp::NumericVector::create(output_years);
  r_aids_deaths_no_art.attr("dim") = Rcpp::NumericVector::create(
      ss.disease_stages, ss.age_groups_hiv, ss.num_genders, output_years);
  r_infections.attr("dim") = Rcpp::NumericVector::create(ss.age_groups_pop, ss.num_genders, output_years);
  r_aids_deaths_art.attr("dim") = Rcpp::NumericVector::create(
      ss.treatment_stages, ss.disease_stages, ss.age_groups_hiv, ss.num_genders, output_years);
  r_art_initiation.attr("dim") = Rcpp::NumericVector::create(
      ss.disease_stages, ss.age_groups_hiv, ss.num_genders, output_years);
  r_hiv_deaths.attr("dim") = Rcpp::NumericVector::create(ss.age_groups_pop, ss.num_genders, output_years);

  std::copy_n(state.total_population.data(), state.total_population.size(),
              REAL(r_total_population));
  std::copy_n(state.natural_deaths.data(), state.natural_deaths.size(),
              REAL(r_natural_deaths));
  std::copy_n(state.hiv_population.data(), state.hiv_population.size(),
              REAL(r_hiv_population));
  std::copy_n(state.hiv_natural_deaths.data(), state.hiv_natural_deaths.size(),
              REAL(r_hiv_natural_deaths));
  std::copy_n(state.hiv_strat_adult.data(), state.hiv_strat_adult.size(),
              REAL(r_hiv_strat_adult));
  std::copy_n(state.art_strat_adult.data(), state.art_strat_adult.size(),
              REAL(r_art_strat_adult));
  std::copy_n(state.births.data(), state.births.size(),
              REAL(r_births));
  std::copy_n(state.aids_deaths_no_art.data(), state.aids_deaths_no_art.size(),
              REAL(r_aids_deaths_no_art));
  std::copy_n(state.infections.data(), state.infections.size(),
              REAL(r_infections));
  std::copy_n(state.aids_deaths_art.data(), state.aids_deaths_art.size(),
              REAL(r_aids_deaths_art));
  std::copy_n(state.art_initiation.data(), state.art_initiation.size(),
              REAL(r_art_initiation));
  std::copy_n(state.hiv_deaths.data(), state.hiv_deaths.size(),
              REAL(r_hiv_deaths));

  return Rcpp::List::create(Rcpp::_["total_population"] = r_total_population,
                            Rcpp::_["births"] = r_births,
                            Rcpp::_["natural_deaths"] = r_natural_deaths,
                            Rcpp::_["hiv_population"] = r_hiv_population,
                            Rcpp::_["hiv_natural_deaths"] = r_hiv_natural_deaths,
                            Rcpp::_["hiv_strat_adult"] = r_hiv_strat_adult,
                            Rcpp::_["art_strat_adult"] = r_art_strat_adult,
                            Rcpp::_["aids_deaths_no_art"] = r_aids_deaths_no_art,
                            Rcpp::_["infections"] = r_infections,
                            Rcpp::_["aids_deaths_art"] = r_aids_deaths_art,
                            Rcpp::_["art_initiation"] = r_art_initiation,
                            Rcpp::_["hiv_deaths"] = r_hiv_deaths);
}
