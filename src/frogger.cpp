#include <Rcpp.h>

#include "frogger.hpp"

// [[Rcpp::export]]
Rcpp::List run_base_model(
  const Rcpp::List data,
  const std::string model_variant,
  const int time_steps,
  const int hiv_steps,
  const std::vector<int> save_steps,
  const bool is_midyear_projection,
  const int t_ART_start
) {
  // TODO write docs, save_steps now must be 0 index based, nothing is nullable, projection_period is "calendar" or "midyear", t_ART_start is 0 based (so substract 1 from R value)
  if (model_variant == "DemographicProjection") {
    return leapfrog::Leapfrog<double, leapfrog::DemographicProjection>::simulate_model(data, time_steps, hiv_steps, save_steps, is_midyear_projection, t_ART_start);
  } else if (model_variant == "HivFullAgeStratification") {
    return leapfrog::Leapfrog<double, leapfrog::HivFullAgeStratification>::simulate_model(data, time_steps, hiv_steps, save_steps, is_midyear_projection, t_ART_start);
  } else if (model_variant == "HivCoarseAgeStratification") {
    return leapfrog::Leapfrog<double, leapfrog::HivCoarseAgeStratification>::simulate_model(data, time_steps, hiv_steps, save_steps, is_midyear_projection, t_ART_start);
  } else if (model_variant == "ChildModel") {
    return leapfrog::Leapfrog<double, leapfrog::ChildModel>::simulate_model(data, time_steps, hiv_steps, save_steps, is_midyear_projection, t_ART_start);
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
