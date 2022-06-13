#include <Rcpp.h>
#include <string>

#include "model_runner.hpp"

//' Simulate base leapfrog model
//'
//' @param demp list of demographic input parameters (TODO: document)
//' @param model_type The type of model, possible options:
//'   * base - Run the base demographic model
//'
//' @details
//' The first year of `survival`, `asfr`, `srb`, and `netmig` is not used. This
//' is assumed ' to apply to the base year population (consistent with
//' Spectrum).
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List demp, const std::string model_type) {
  Rcpp::NumericVector survival = demp["survival"];
  Rcpp::Dimension d = survival.attr("dim");
  const size_t proj_years =
      2;  // TODO: Fix temporarily for testing should be d[2]
  const int num_genders = 2;
  const int age_groups_pop = 81;
  const int fertility_first_age_group = 15;
  const int age_groups_fert = 35;

  // allocate memory for return object
  Rcpp::NumericVector total_population(age_groups_pop * num_genders);
  total_population.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);
  double births;

  Rcpp::NumericVector natural_deaths(age_groups_pop * num_genders);
  natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);

  leapfrog_sim<double>(
      REAL(demp["basepop"]), REAL(demp["Sx"]), REAL(demp["netmigr_adj"]),
      REAL(demp["asfr"]), REAL(demp["births_sex_prop"]), model_type,
      num_genders, age_groups_pop, fertility_first_age_group, age_groups_fert,
      proj_years, REAL(total_population), &births, REAL(natural_deaths));

  Rcpp::List ret = Rcpp::List::create(
      Rcpp::_["total_population"] = total_population,
      Rcpp::_["births"] = births, Rcpp::_["natural_deaths"] = natural_deaths);

  return ret;
}