#include <Rcpp.h>

#include "code.hpp"
#include "frogger.hpp"

int get_simulation_years(const Rcpp::List demp, SEXP r_sim_years) {
  Rcpp::NumericVector Sx = demp["Sx"];
  Rcpp::Dimension d = Sx.attr("dim");
  const int max_sim_years = d[2];
  if (r_sim_years == R_NilValue) {
    return max_sim_years;
  }
  auto sim_years = INTEGER(r_sim_years)[0];
  if (sim_years > max_sim_years) {
    Rcpp::stop("Too long");
  }
  return sim_years;
}

// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List demp, SEXP sim_years) {
  const int proj_years = get_simulation_years(demp, sim_years);
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

  auto state = model_runner(proj_years, params);

  Rcpp::NumericVector r_total_population(age_groups_pop * num_genders);
  Rcpp::NumericVector r_births(1);
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
