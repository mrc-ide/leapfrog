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
    "ChildModel",
    "GBD"
  };
}

template<typename ModelVariant>
using LeapfrogR = leapfrog::Leapfrog<leapfrog::R, double, ModelVariant>;

auto get_opts_r(
  const Rcpp::List parameters,
  const std::vector<int> output_years
) {
  const int t_art_start = Rcpp::as<int>(parameters["t_ART_start"]);
  const int hts_per_year = Rcpp::as<int>(parameters["hts_per_year"]);
  const int proj_start_year = Rcpp::as<int>(
    parameters["projection_start_year"]);
  const std::string projection_period = Rcpp::as<std::string>(
    parameters["projection_period"]);

  return leapfrog::get_opts<double>(
    hts_per_year, t_art_start, projection_period, proj_start_year, output_years
  );
}



template<typename ModelVariant>
auto build_output_r(
  typename LeapfrogR<ModelVariant>::OutputState output_state,
  const std::vector<int> output_years
) {
  using LF = LeapfrogR<ModelVariant>;

  const int output_size = LF::Cfg::get_build_output_size(0);
  Rcpp::List ret(output_size);
  Rcpp::CharacterVector names(output_size);
  LF::Cfg::build_output(0, output_state, ret, names, output_years.size());
  ret.attr("names") = names;

  return ret;
}

template<typename ModelVariant>
auto build_output_r(
  typename LeapfrogR<ModelVariant>::State state
) {
  using LF = LeapfrogR<ModelVariant>;

  const int output_size = LF::Cfg::get_build_output_size(0);
  Rcpp::List ret(output_size);
  Rcpp::CharacterVector names(output_size);
  LF::Cfg::build_output_single_year(0, state, ret, names);
  ret.attr("names") = names;

  return ret;
}



template<typename ModelVariant>
Rcpp::List simulate_model(
  const Rcpp::List parameters,
  const std::vector<int> output_years
) {
  using LF = LeapfrogR<ModelVariant>;

  const auto opts = get_opts_r(parameters, output_years);
  const auto pars = LF::Cfg::get_pars(parameters, opts);

  auto state = LF::run_model(pars, opts, output_years);

  return build_output_r<ModelVariant>(state, output_years);
}

template<typename ModelVariant>
Rcpp::List simulate_model(
  const Rcpp::List parameters,
  const std::vector<int> output_years,
  Rcpp::List initial_state_data,
  int simulation_start_year
) {
  using LF = LeapfrogR<ModelVariant>;

  const auto opts = get_opts_r(parameters, output_years);
  const auto pars = LF::Cfg::get_pars(parameters, opts);

  typename LF::State initial_state = LF::Cfg::get_initial_state(initial_state_data);
  auto state = LF::run_model_from_state(pars, opts, initial_state, simulation_start_year, output_years);

  return build_output_r<ModelVariant>(state, output_years);
}

template<typename ModelVariant>
Rcpp::List simulate_model(
  const Rcpp::List parameters,
  Rcpp::List initial_state_data,
  int simulation_start_year
) {
  using LF = LeapfrogR<ModelVariant>;

  const auto opts = get_opts_r(parameters, { simulation_start_year + 1 });
  const auto pars = LF::Cfg::get_pars(parameters, opts);

  typename LF::State initial_state = LF::Cfg::get_initial_state(initial_state_data);
  auto state = LF::run_model_single_year(pars, opts, initial_state, simulation_start_year);

  return build_output_r<ModelVariant>(state);
}



template<typename ...Args>
auto sim_model(const std::string configuration, Args&&... args) {
  if (configuration == "DemographicProjection") {
    return simulate_model<leapfrog::DemographicProjection>(std::forward<Args>(args)...);
  } else if (configuration == "HivFullAgeStratification") {
    return simulate_model<leapfrog::HivFullAgeStratification>(std::forward<Args>(args)...);
  } else if (configuration == "HivCoarseAgeStratification") {
    return simulate_model<leapfrog::HivCoarseAgeStratification>(std::forward<Args>(args)...);
  } else if (configuration == "ChildModel") {
    return simulate_model<leapfrog::ChildModel>(std::forward<Args>(args)...);
  } else if (configuration == "GBD") {
    return simulate_model<leapfrog::GBD>(std::forward<Args>(args)...);    
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
  return sim_model(configuration, parameters, output_years);
}

// [[Rcpp::export]]
Rcpp::List run_base_model_from_state(
  const Rcpp::List parameters,
  const std::string configuration,
  const Rcpp::List initial_state,
  int simulation_start_year,
  const std::vector<int> output_years
) {
  return sim_model(configuration, parameters, output_years, initial_state, simulation_start_year);
}

// [[Rcpp::export]]
Rcpp::List run_base_model_single_year(
  const Rcpp::List parameters,
  const std::string configuration,
  const Rcpp::List initial_state,
  int simulation_start_year
) {
  return sim_model(configuration, parameters, initial_state, simulation_start_year);
}
