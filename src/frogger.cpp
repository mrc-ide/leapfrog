#include <Rcpp.h>

#include "frogger.hpp"
#include "intermediate_data.hpp"
#include "model_variants.hpp"
#include "state_space.hpp"
#include "generated/model_input.hpp"
#include "generated/model_output.hpp"
#include "r_utils.hpp"

template<typename ModelVariant>
Rcpp::List simulate_model(const leapfrog::StateSpace<ModelVariant> ss,
                          const Rcpp::List data,
                          const int proj_years,
                          const int hiv_steps,
                          const std::vector<int> save_steps) {
  const leapfrog::Options<double> opts = {
      hiv_steps,
      Rcpp::as<int>(data["t_ART_start"]) - 1,
      ss.base.hAG
  };

  const auto params = setup_model_params<ModelVariant, double>(data, opts, proj_years);

  auto state = leapfrog::run_model<ModelVariant, double>(proj_years, save_steps, params);

  auto ret = build_r_output<ModelVariant, double>(state, save_steps);

  return ret;
}

// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List data,
                          const int proj_years,
                          const int hiv_steps,
                          const std::vector<int> save_steps,
                          std::string model_variant) {
  Rcpp::List ret;
  if (model_variant == "ChildModel") {
    constexpr auto ss = leapfrog::StateSpace<leapfrog::ChildModel>();
    ret = simulate_model<leapfrog::ChildModel>(ss, data, proj_years, hiv_steps, save_steps);
  } else if (model_variant == "BaseModelFullAgeStratification") {
    constexpr auto ss = leapfrog::StateSpace<leapfrog::BaseModelFullAgeStratification>();
    ret = simulate_model<leapfrog::BaseModelFullAgeStratification>(ss, data, proj_years, hiv_steps, save_steps);
  } else if (model_variant == "BaseModelCoarseAgeStratification") {
    constexpr auto ss = leapfrog::StateSpace<leapfrog::BaseModelCoarseAgeStratification>();
    ret = simulate_model<leapfrog::BaseModelCoarseAgeStratification>(ss, data, proj_years, hiv_steps, save_steps);
  } else {
    Rcpp::stop("Invalid model variant " + model_variant + " must be one of " +
               "'BaseModelFullAgeStratification', 'BaseModelCoarseAgeStratification' or 'ChildModel'");
  }

  return ret;
}
