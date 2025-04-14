#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>

#include "frogger.hpp"
#include "generated/r_interface/r_adapters.hpp"


//' List the avaialble model configurations
//'
//' @return List of available model configurations
//' @export
// [[Rcpp::export]]
std::vector<std::string> list_model_configurations() {
  return std::vector<std::string>{
    "DemographicProjection",
    "HivFullAgeStratification",
    "HivCoarseAgeStratification",
    "ChildModel"
  };
}

template<typename real_type, typename ModelVariant>
Rcpp::List simulate_model(
  const Rcpp::List parameters,
  const std::vector<int> output_years
) {
  using LF = leapfrog::Leapfrog<leapfrog::R, real_type, ModelVariant>;

  const int t_art_start = Rcpp::as<int>(parameters["t_ART_start"]);
  const int hts_per_year = Rcpp::as<int>(parameters["hts_per_year"]);
  const int proj_start_year = Rcpp::as<int>(
    parameters["projection_start_year"]);
  const std::string projection_period = Rcpp::as<std::string>(
    parameters["projection_period"]);

  const auto opts = leapfrog::get_opts<real_type>(hts_per_year,
                                                  t_art_start,
                                                  projection_period,
                                                  proj_start_year,
                                                  output_years);
  const auto pars = LF::Cfg::get_pars(parameters, opts);

  auto state = LF::run_model(pars, opts, output_years);

  const int output_size = LF::Cfg::get_build_output_size(0);
  Rcpp::List ret(output_size);
  Rcpp::CharacterVector names(output_size);
  LF::Cfg::build_output(0, state, ret, names, output_years.size());
  ret.attr("names") = names;

  return ret;
}

// [[Rcpp::export]]
Rcpp::List run_base_model(
  const Rcpp::List parameters,
  const std::string configuration,
  const std::vector<int> output_years
) {
  if (configuration == "DemographicProjection") {
    return simulate_model<double, leapfrog::DemographicProjection>(parameters, output_years);
  } else if (configuration == "HivFullAgeStratification") {
    return simulate_model<double, leapfrog::HivFullAgeStratification>(parameters, output_years);
  } else if (configuration == "HivCoarseAgeStratification") {
    return simulate_model<double, leapfrog::HivCoarseAgeStratification>(parameters, output_years);
  } else if (configuration == "ChildModel") {
    return simulate_model<double, leapfrog::ChildModel>(parameters, output_years);
  } else {
    const auto available_variants = list_model_configurations();
    std::ostringstream oss;
    oss << "Invalid configuration: '" << configuration
        << "'. It must be one of: ";

    for (size_t i = 0; i < available_variants.size(); ++i) {
      oss << "'" << available_variants[i] << "'";
      if (i != available_variants.size() - 1) {
        oss << ", ";
      } else {
        oss << ".";
      }
    }

    throw std::runtime_error(oss.str());
  }
}
