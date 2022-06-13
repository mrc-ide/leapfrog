#include <Rcpp.h>

#include "code.hpp"

// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List demp) {
  const size_t proj_years =
      2;  // TODO: Fix temporarily for testing should be d[2]
  const int num_genders = 2;
  const int age_groups_pop = 81;
  const int fertility_first_age_group = 15;
  const int age_groups_fert = 35;

  TensorMapX2T<double> base_pop(REAL(demp["basepop"]), age_groups_pop,
                                num_genders);
  TensorMapX2T<double> survival(REAL(demp["Sx"]), age_groups_pop + 1,
                                num_genders);
  TensorMapX2T<double> net_migration(REAL(demp["netmigr_adj"]), age_groups_pop,
                                     num_genders);
  TensorMapX1T<double> age_sex_fertility_ratio(REAL(demp["asfr"]),
                                               age_groups_fert);
  TensorMapX1T<double> births_sex_prop(REAL(demp["births_sex_prop"]),
                                       num_genders);

  Parameters<double> params = {num_genders,
                               age_groups_pop,
                               fertility_first_age_group,
                               age_groups_fert,
                               base_pop,
                               survival,
                               net_migration,
                               age_sex_fertility_ratio,
                               births_sex_prop};

  State<double> state(age_groups_pop, num_genders);

  // leapfrog_sim<double>(
  //     REAL(demp["basepop"]), REAL(demp["Sx"]), REAL(demp["netmigr_adj"]),
  //     REAL(demp["asfr"]), REAL(demp["births_sex_prop"]), num_genders,
  //     age_groups_pop, fertility_first_age_group, age_groups_fert,
  //     proj_years, REAL(total_population), REAL(births),
  //     REAL(natural_deaths));

  Rcpp::NumericVector r_total_population(age_groups_pop * num_genders);
  Rcpp::NumericVector r_births;
  Rcpp::NumericVector r_natural_deaths(age_groups_pop * num_genders);
  r_total_population.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);
  r_natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);

  std::copy_n(state.total_population.data(), state.total_population.size(),
              REAL(r_total_population));
  std::copy_n(state.natural_deaths.data(), state.natural_deaths.size(),
              REAL(r_natural_deaths));
  REAL(r_births)[0] = state.births;

  Rcpp::List ret =
      Rcpp::List::create(Rcpp::_["total_population"] = r_total_population,
                         Rcpp::_["births"] = r_births,
                         Rcpp::_["natural_deaths"] = r_natural_deaths);
  return ret;
}
