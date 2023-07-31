#include <Rcpp.h>

#include "frogger.hpp"
#include "types.hpp"
#include "state_space.hpp"
#include "model_setup.h"

int transform_simulation_years(const Rcpp::List demp, SEXP r_sim_years) {
  Rcpp::NumericVector Sx = demp["Sx"];
  Rcpp::Dimension d = Sx.attr("dim");
  // Simulation initialises state from first years input data (index 0)
  // then runs for each year simulating this years (i) data using previous years
  // state (i - 1) and this years input data (i). So -1 off index for max years
  // to simulate as index 0 used for initial state
  const int max_sim_years = d[2] - 1;
  if (r_sim_years == R_NilValue) {
    return max_sim_years;
  }
  auto sim_years = INTEGER(r_sim_years)[0];
  if (sim_years > max_sim_years) {
    Rcpp::stop("No of years > max years of " + std::to_string(max_sim_years));
  }
  return sim_years;
}

int transform_hiv_steps_per_year(SEXP r_hiv_steps_per_year) {
  int hiv_steps_per_year;
  if (r_hiv_steps_per_year == R_NilValue) {
    hiv_steps_per_year = 10;
  } else {
    hiv_steps_per_year = INTEGER(r_hiv_steps_per_year)[0];
  }
  return hiv_steps_per_year;
}

template <leapfrog::HivAgeStratification S>
leapfrog::StateSpace<S> get_state_space(const std::string hiv_age_stratification) {
  leapfrog::StateSpace<S> ss;
  if (hiv_age_stratification == "full") {
    ss = leapfrog::StateSpace<leapfrog::HivAgeStratification::full>();
  } else {
    // We've already validated that hiv_age_stratification is one of "full" or "coarse"
    ss = leapfrog::StateSpace<leapfrog::HivAgeStratification::coarse>();
  }
  return ss;
}

std::vector<int> transform_output_steps(Rcpp::NumericVector output_steps) {
  return Rcpp::as<std::vector<int>>(output_steps);
}

template <leapfrog::HivAgeStratification S>
Rcpp::List fit_model(const leapfrog::StateSpace<S> ss,
                     const Rcpp::List data,
                     const int proj_years,
                     const int hiv_steps,
                     const bool run_child_model,
                     const std::vector<int> save_steps) {
  const leapfrog::Options<double> opts = {
      hiv_steps,
      Rcpp::as<int>(data["t_ART_start"]) - 1,
      ss.age_groups_hiv,
      Rcpp::as<int>(data["scale_cd4_mort"]),
      Rcpp::as<double>(data["art_alloc_mxweight"]),
      run_child_model
  };

  const leapfrog::Parameters<double> params = setup_model_params<double, S>(data, opts, proj_years);

  auto state = leapfrog::run_model<S, double>(proj_years, save_steps, params);

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
  Rcpp::NumericVector r_hc1_hiv_pop(ss.hc1_disease_stages * ss.hTM * ss.hc1_age_groups * ss.num_genders * output_years);
  Rcpp::NumericVector r_hc2_hiv_pop(ss.hc2_disease_stages * ss.hTM * ss.hc2_age_groups * ss.num_genders * output_years);
  Rcpp::NumericVector r_hc1_art_pop(ss.treatment_stages * ss.hc1_disease_stages * ss.hc1_age_groups * ss.num_genders * output_years);
  Rcpp::NumericVector r_hc2_art_pop(ss.treatment_stages * ss.hc2_disease_stages * ss.hc2_age_groups * ss.num_genders * output_years);
  Rcpp::NumericVector r_hc1_art_aids_deaths(ss.treatment_stages * ss.hc1_disease_stages * ss.hc1_age_groups * ss.num_genders * output_years);
  Rcpp::NumericVector r_hc2_art_aids_deaths(ss.treatment_stages * ss.hc2_disease_stages * ss.hc2_age_groups * ss.num_genders * output_years);
  Rcpp::NumericVector r_hc1_noart_aids_deaths(ss.hc1_disease_stages * ss.hTM * ss.hc1_age_groups * ss.num_genders * output_years);
  Rcpp::NumericVector r_hc2_noart_aids_deaths(ss.hc2_disease_stages * ss.hTM * ss.hc2_age_groups * ss.num_genders * output_years);

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
  r_hc1_hiv_pop.attr("dim") = Rcpp::NumericVector::create(ss.hc1_disease_stages, ss.hTM, ss.hc1_age_groups, ss.num_genders,
                                                         output_years);
  r_hc2_hiv_pop.attr("dim") = Rcpp::NumericVector::create(ss.hc2_disease_stages, ss.hTM, ss.hc2_age_groups, ss.num_genders,
                    output_years);
  r_hc1_art_pop.attr("dim") = Rcpp::NumericVector::create(ss.treatment_stages, ss.hc1_disease_stages, ss.hc1_age_groups, ss.num_genders,
                     output_years);
  r_hc2_art_pop.attr("dim") = Rcpp::NumericVector::create(ss.treatment_stages, ss.hc2_disease_stages, ss.hc2_age_groups, ss.num_genders,
                     output_years);
  r_hc1_art_aids_deaths.attr("dim") = Rcpp::NumericVector::create(ss.treatment_stages, ss.hc1_disease_stages, ss.hc1_age_groups, ss.num_genders,
                     output_years);
  r_hc2_art_aids_deaths.attr("dim") = Rcpp::NumericVector::create(ss.treatment_stages, ss.hc2_disease_stages, ss.hc2_age_groups, ss.num_genders,
                     output_years);
  r_hc1_noart_aids_deaths.attr("dim") = Rcpp::NumericVector::create( ss.hc1_disease_stages, ss.hTM, ss.hc1_age_groups, ss.num_genders,
                             output_years);
  r_hc2_noart_aids_deaths.attr("dim") = Rcpp::NumericVector::create(ss.hc2_disease_stages, ss.hTM, ss.hc2_age_groups, ss.num_genders,
                             output_years);

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
  std::copy_n(state.hc1_hiv_pop.data(), state.hc1_hiv_pop.size(),
              REAL(r_hc1_hiv_pop));
  std::copy_n(state.hc2_hiv_pop.data(), state.hc2_hiv_pop.size(),
              REAL(r_hc2_hiv_pop));
  std::copy_n(state.hc1_art_pop.data(), state.hc1_art_pop.size(),
              REAL(r_hc1_art_pop));
  std::copy_n(state.hc2_art_pop.data(), state.hc2_art_pop.size(),
              REAL(r_hc2_art_pop));
  std::copy_n(state.hc1_art_aids_deaths.data(), state.hc1_art_aids_deaths.size(),
              REAL(r_hc1_art_aids_deaths));
  std::copy_n(state.hc2_art_aids_deaths.data(), state.hc2_art_aids_deaths.size(),
              REAL(r_hc2_art_aids_deaths));
  std::copy_n(state.hc1_noart_aids_deaths.data(), state.hc1_noart_aids_deaths.size(),
              REAL(r_hc1_noart_aids_deaths));
  std::copy_n(state.hc2_noart_aids_deaths.data(), state.hc2_noart_aids_deaths.size(),
              REAL(r_hc2_noart_aids_deaths));

  Rcpp::List ret =
      Rcpp::List::create(Rcpp::_["total_population"] = r_total_population,
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
                         Rcpp::_["hiv_deaths"] = r_hiv_deaths,
                         Rcpp::_["hc1_hiv_pop"] = r_hc1_hiv_pop,
                         Rcpp::_["hc2_hiv_pop"] = r_hc2_hiv_pop,
                         Rcpp::_["hc1_art_pop"] = r_hc1_art_pop,
                         Rcpp::_["hc2_art_pop"] = r_hc2_art_pop,
                         Rcpp::_["hc1_art_aids_deaths"] = r_hc1_art_aids_deaths,
                         Rcpp::_["hc2_art_aids_deaths"] = r_hc2_art_aids_deaths,
                         Rcpp::_["hc1_noart_aids_deaths"] = r_hc1_noart_aids_deaths,
                         Rcpp::_["hc2_noart_aids_deaths"] = r_hc2_noart_aids_deaths);
  return ret;
}

// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List data,
                          SEXP sim_years,
                          SEXP hiv_steps_per_year,
                          Rcpp::NumericVector output_steps,
                          std::string hiv_age_stratification,
                          bool run_child_model) {
  const int proj_years = transform_simulation_years(data, sim_years);
  const std::vector<int> save_steps = transform_output_steps(output_steps);
  const int hiv_steps = transform_hiv_steps_per_year(hiv_steps_per_year);


  Rcpp::List ret;
  if (hiv_age_stratification == "full") {
    constexpr auto ss = leapfrog::StateSpace<leapfrog::HivAgeStratification::full>();
    ret = fit_model<leapfrog::HivAgeStratification::full>(ss, data,
                                                          proj_years, hiv_steps, run_child_model, save_steps);
  } else {
    // We've already validated that hiv_age_stratification is one of "full" or "coarse"
    constexpr auto ss = leapfrog::StateSpace<leapfrog::HivAgeStratification::coarse>();
    ret = fit_model<leapfrog::HivAgeStratification::coarse>(ss, data,
                                                            proj_years, hiv_steps, run_child_model, save_steps);
  }

  return ret;
}
