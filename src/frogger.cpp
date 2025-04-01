#include <Rcpp.h>

#include "frogger.hpp"
#include "generated/config_mixer.hpp"
#include "generated/concepts.hpp"
#include "utils/language_types.hpp"

template <typename real_type, typename ModelVariant>
Rcpp::List run_frogger(
  const Rcpp::List data,
  const std::string model_variant,
  const int time_steps,
  const int hiv_steps,
  const std::vector<int> save_steps,
  const bool is_midyear_projection,
  const int t_ART_start
) {
  const int output_size = leapfrog::ConfigMixed<real_type, ModelVariant>::get_build_output_size(0);
  Rcpp::List out_data(output_size);
  Rcpp::CharacterVector names(output_size);
  out_data.attr("names") = names;
  OutputData ret = {
    out_data,
    names
  };
  leapfrog::Leapfrog<real_type, ModelVariant>::simulate_model(data, time_steps, hiv_steps, save_steps, is_midyear_projection, t_ART_start, ret);
  return ret.data;
}

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
  // TODO Mantra write docs, save_steps now must be 0 index based, nothing is nullable, projection_period is "calendar" or "midyear", t_ART_start is 0 based (so substract 1 from R value)
  if (model_variant == "DemographicProjection") {
    return run_frogger<double, leapfrog::DemographicProjection>(data, model_variant, time_steps, hiv_steps, save_steps, is_midyear_projection, t_ART_start);
  } else if (model_variant == "HivFullAgeStratification") {
    return run_frogger<double, leapfrog::HivFullAgeStratification>(data, model_variant, time_steps, hiv_steps, save_steps, is_midyear_projection, t_ART_start);
  } else if (model_variant == "HivCoarseAgeStratification") {
    return run_frogger<double, leapfrog::HivCoarseAgeStratification>(data, model_variant, time_steps, hiv_steps, save_steps, is_midyear_projection, t_ART_start);
  } else if (model_variant == "ChildModel") {
    return run_frogger<double, leapfrog::ChildModel>(data, model_variant, time_steps, hiv_steps, save_steps, is_midyear_projection, t_ART_start);
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
