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

template<typename ModelVariant>
Rcpp::List simulate_model(
  const Rcpp::List parameters,
  const std::vector<int> output_years,
  const std::optional<Rcpp::List> initial_state_data = std::nullopt,
  std::optional<int> start_from_year_nullable = std::nullopt
) {
  using LF = leapfrog::Leapfrog<leapfrog::R, double, ModelVariant>;

  const int t_art_start = Rcpp::as<int>(parameters["t_ART_start"]);
  const bool is_midyear_projection = Rcpp::as<bool>(
    parameters["is_midyear_projection"]);
  int hts_per_year = 10;
  if (parameters.containsElementNamed("hts_per_year")) {
    hts_per_year = Rcpp::as<int>(parameters["hts_per_year"]);
  }
  const int proj_start_year = Rcpp::as<int>(
    parameters["projection_start_year"]);

  const auto opts = leapfrog::get_opts<double>(
    hts_per_year, t_art_start, is_midyear_projection, proj_start_year, output_years
  );

  const auto pars = LF::Cfg::get_pars(parameters, opts);
  
  typename LF::State initial_state = {};
  if (initial_state_data) {
    initial_state = LF::Cfg::get_initial_state(initial_state_data.value());
  } else {
    initial_state.reset();
    if constexpr (ModelVariant::run_demographic_projection) {
      initial_state.dp.p_total_pop = pars.dp.base_pop;
    }
  }

  int start_from_year;
  if (start_from_year_nullable) {
    start_from_year = start_from_year_nullable.value();
  } else {
    start_from_year = opts.proj_start_year;
  }

  auto state = LF::run_model(pars, opts, initial_state, start_from_year, output_years);

  const int output_size = LF::Cfg::get_build_output_size(0);
  Rcpp::List ret(output_size);
  Rcpp::CharacterVector names(output_size);
  LF::Cfg::build_output(0, state, ret, names, output_years.size());
  ret.attr("names") = names;

  return ret;
}

auto get_sim_model(const std::string configuration) {
  if (configuration == "DemographicProjection") {
    return simulate_model<leapfrog::DemographicProjection>;
  } else if (configuration == "HivFullAgeStratification") {
    return simulate_model<leapfrog::HivFullAgeStratification>;
  } else if (configuration == "HivCoarseAgeStratification") {
    return simulate_model<leapfrog::HivCoarseAgeStratification>;
  } else if (configuration == "ChildModel") {
    return simulate_model<leapfrog::ChildModel>;
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
};

// [[Rcpp::export]]
Rcpp::List run_base_model(
  const Rcpp::List parameters,
  const std::string configuration,
  const std::vector<int> output_years
) {
  return get_sim_model(configuration)(parameters, output_years, std::nullopt, std::nullopt);
}

// [[Rcpp::export]]
Rcpp::List run_base_model_with_initial_state(
  const Rcpp::List parameters,
  const std::string configuration,
  const std::vector<int> output_years,
  const Rcpp::List initial_state,
  int start_from_year
) {
  return get_sim_model(configuration)(parameters, output_years, initial_state, start_from_year);
}
