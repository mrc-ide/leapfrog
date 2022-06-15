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

// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List demp, SEXP sim_years) {
  const int proj_years = get_simulation_years(demp, sim_years);
  const int num_genders = 2;
  const int age_groups_pop = 81;
  const int fertility_first_age_group = 15;
  const int age_groups_fert = 35;

  TensorMap2<double> base_pop(REAL(demp["basepop"]), age_groups_pop,
                              num_genders);
  // Survival has size age_groups_pop + 1 as this is the probability of
  // surviving between ages, so from 0 to 1, 1 to 2, ..., 79 to 80+ and
  // 80+ to 80+
  TensorMap3<double> survival(REAL(demp["Sx"]), age_groups_pop + 1, num_genders,
                              proj_years);
  TensorMap3<double> net_migration(REAL(demp["netmigr_adj"]), age_groups_pop,
                                   num_genders, proj_years);
  TensorMap2<double> age_sex_fertility_ratio(REAL(demp["asfr"]),
                                             age_groups_fert, proj_years);
  TensorMap2<double> births_sex_prop(REAL(demp["births_sex_prop"]), num_genders,
                                     proj_years);

  Parameters<double> params = {num_genders,
                               age_groups_pop,
                               fertility_first_age_group,
                               age_groups_fert,
                               base_pop,
                               survival,
                               net_migration,
                               age_sex_fertility_ratio,
                               births_sex_prop};

  auto state = run_model(proj_years, params);

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
