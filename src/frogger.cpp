#include <Rcpp.h>

#include "frogger.hpp"
#include "types.hpp"

int get_simulation_years(const Rcpp::List demp, SEXP r_sim_years) {
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

TensorMap1<int> get_age_groups_hiv_span(const Rcpp::List projection_parameters,
                                        std::string hiv_age_stratification) {
  int age_groups_hiv;
  SEXP data;
  if (hiv_age_stratification == "full") {
    age_groups_hiv = 66;
    data = projection_parameters["hAG_SPAN_full"];
  } else if (hiv_age_stratification == "coarse") {
    age_groups_hiv = 9;
    data = projection_parameters["hAG_SPAN_coarse"];
  } else {
    Rcpp::stop(
        "Invalid HIV age stratification must be 'full' or 'coarse' got '%s'.",
        hiv_age_stratification);
  }
  if (LENGTH(data) != age_groups_hiv) {
    Rcpp::stop("Invalid size of data, expected %d got %d", age_groups_hiv,
               LENGTH(data));
  }
  TensorMap1<int> age_groups_hiv_span(INTEGER(data), age_groups_hiv);
  return age_groups_hiv_span;
}

// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List data,
                          const Rcpp::List projection_parameters,
                          SEXP sim_years,
                          std::string hiv_age_stratification = "full") {
  const int proj_years = get_simulation_years(data, sim_years);
  const int num_genders = 2;
  const int age_groups_pop = 81;
  const int fertility_first_age_group = 15;
  const int age_groups_fert = 35;
  const int disease_stages = 7;
  const int treatment_stages = 3;
  const int hiv_adult_first_age_group = 15;
  const int time_art_start = *INTEGER(projection_parameters["t_ART_start"]) -
                             1;  // 0-based indexing vs R 1-based
  auto age_groups_hiv_span =
      get_age_groups_hiv_span(projection_parameters, hiv_age_stratification);
  int age_groups_hiv = static_cast<int>(age_groups_hiv_span.size());

  TensorMap2<double> base_pop(REAL(data["basepop"]), age_groups_pop,
                              num_genders);
  // Survival has size age_groups_pop + 1 as this is the probability of
  // surviving between ages, so from 0 to 1, 1 to 2, ..., 79 to 80+ and
  // 80+ to 80+
  TensorMap3<double> survival(REAL(data["Sx"]), age_groups_pop + 1, num_genders,
                              proj_years);
  TensorMap3<double> net_migration(REAL(data["netmigr_adj"]), age_groups_pop,
                                   num_genders, proj_years);
  TensorMap2<double> age_sex_fertility_ratio(REAL(data["asfr"]),
                                             age_groups_fert, proj_years);
  TensorMap2<double> births_sex_prop(REAL(data["births_sex_prop"]), num_genders,
                                     proj_years);

  Parameters<double> params = {num_genders,
                               age_groups_pop,
                               fertility_first_age_group,
                               age_groups_fert,
                               age_groups_hiv,
                               disease_stages,
                               hiv_adult_first_age_group,
                               treatment_stages,
                               time_art_start,
                               age_groups_hiv_span,
                               base_pop,
                               survival,
                               net_migration,
                               age_sex_fertility_ratio,
                               births_sex_prop};

  auto state = run_model(proj_years, params);

  Rcpp::NumericVector r_total_population(age_groups_pop * num_genders);
  Rcpp::NumericVector r_births(1);
  Rcpp::NumericVector r_natural_deaths(age_groups_pop * num_genders);
  Rcpp::NumericVector r_hiv_population(age_groups_pop * num_genders);
  Rcpp::NumericVector r_hiv_natural_deaths(age_groups_pop * num_genders);
  Rcpp::NumericVector r_hiv_strat_adult(disease_stages * age_groups_hiv *
                                        num_genders);
  Rcpp::NumericVector r_art_strat_adult(treatment_stages * disease_stages *
                                        age_groups_hiv * num_genders);
  r_total_population.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);
  r_natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);
  r_hiv_population.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);
  r_hiv_natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);
  r_hiv_strat_adult.attr("dim") =
      Rcpp::NumericVector::create(disease_stages, age_groups_hiv, num_genders);
  r_art_strat_adult.attr("dim") = Rcpp::NumericVector::create(
      treatment_stages, disease_stages, age_groups_hiv, num_genders);

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
  REAL(r_births)[0] = state.births;

  Rcpp::List ret =
      Rcpp::List::create(Rcpp::_["total_population"] = r_total_population,
                         Rcpp::_["births"] = r_births,
                         Rcpp::_["natural_deaths"] = r_natural_deaths,
                         Rcpp::_["hiv_population"] = r_hiv_population,
                         Rcpp::_["hiv_natural_deaths"] = r_hiv_natural_deaths,
                         Rcpp::_["hiv_strat_adult"] = r_hiv_strat_adult,
                         Rcpp::_["art_strat_adult"] = r_art_strat_adult);
  return ret;
}
