#include <Rcpp.h>

#include "frogger.hpp"
#include "generated/r_interface/r_adapters.hpp"

template<typename real_type, typename ModelVariant>
Rcpp::List simulate_model(
  const Rcpp::List parameters,
  const std::vector<int> output_years
) {
  using LF = leapfrog::Leapfrog<leapfrog::R, real_type, ModelVariant>;

  const int t_art_start = Rcpp::as<int>(parameters["t_ART_start"]);
  const bool is_midyear_projection = Rcpp::as<bool>(
    parameters["is_midyear_projection"]);
  int hts_per_year = 10;
  if (parameters.containsElementNamed("hts_per_year")) {
    hts_per_year = Rcpp::as<int>(parameters["hts_per_year"]);
  }
  const int proj_start_year = Rcpp::as<int>(
    parameters["projection_start_year"]);

  const auto opts = leapfrog::get_opts<real_type>(hts_per_year,
                                                  t_art_start,
                                                  is_midyear_projection,
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
  const std::string model_variant,
  const std::vector<int> output_years
) {
  // TODO Mantra write docs, save_steps now must be 0 index based, nothing is nullable, projection_period is "calendar" or "midyear", t_ART_start is 0 based (so substract 1 from R value)
  if (model_variant == "DemographicProjection") {
    return simulate_model<double, leapfrog::DemographicProjection>(parameters, output_years);
  } else if (model_variant == "HivFullAgeStratification") {
    return simulate_model<double, leapfrog::HivFullAgeStratification>(parameters, output_years);
  } else if (model_variant == "HivCoarseAgeStratification") {
    return simulate_model<double, leapfrog::HivCoarseAgeStratification>(parameters, output_years);
  } else if (model_variant == "ChildModel") {
    return simulate_model<double, leapfrog::ChildModel>(parameters, output_years);
  } else {
    throw std::runtime_error(
      "Invalid model_variant: " + model_variant +
      ". It must be one of" +
      "'DemographicProjection', " +
      "'HivFullAgeStratification', " +
      "'HivCoarseAgeStratification', " +
      "'ChildModel'"
    );
  }
}
