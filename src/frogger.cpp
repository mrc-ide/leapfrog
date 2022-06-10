#include <Rcpp.h>
#include <string>

#include "model_runner.hpp"

//' Add two numbers.
//'
//' @param a First number.
//' @param b Second number.
//' @return The sum of two values.
//' @export
// [[Rcpp::export]]
int adder(int a, int b) {
  return add(a, b);
}

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
Rcpp::List run_base_model(const Rcpp::List& demp,
                          const std::string model_type) {
  Rcpp::NumericVector survival = demp["survival"];
  Rcpp::Dimension d = survival.attr("dim");
  const size_t proj_years =
      2;  // TODO: Fix temporarily for testing should be d[2]
  const int NUM_GENDERS = 2;
  const int AGE_GROUPS_POP = 81;
  const int FERTILITY_FIRST_AGE_GROUP = 15;
  const int AGE_GROUPS_FERT = 35;

  // allocate memory for return object
  Rcpp::NumericVector total_population(AGE_GROUPS_POP * NUM_GENDERS *
                                       proj_years);
  total_population.attr("dim") =
      Rcpp::NumericVector::create(AGE_GROUPS_POP, NUM_GENDERS, proj_years);
  Rcpp::NumericVector births(proj_years);

  Rcpp::NumericVector natural_deaths(AGE_GROUPS_POP * NUM_GENDERS * proj_years);
  natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(AGE_GROUPS_POP, NUM_GENDERS, proj_years);

  leapfrog_sim<double, NUM_GENDERS, AGE_GROUPS_POP, FERTILITY_FIRST_AGE_GROUP,
               AGE_GROUPS_FERT>(
      REAL(demp["basepop"]), REAL(demp["Sx"]), REAL(demp["netmigr_adj"]),
      REAL(demp["asfr"]), REAL(demp["births_sex_prop"]), model_type, proj_years,
      REAL(total_population), REAL(births), REAL(natural_deaths));

  Rcpp::List ret = Rcpp::List::create(
      Rcpp::_["total_population"] = total_population,
      Rcpp::_["births"] = births, Rcpp::_["natural_deaths"] = natural_deaths);

  return ret;
}