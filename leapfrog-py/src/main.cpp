#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <vector>
#include <sstream>

#include "../../r-package/inst/include/frogger.hpp"
#include "../../r-package/inst/include/generated/py_interface/py_adapters.hpp"

template<typename ModelVariant>
using LeapfrogPy = leapfrog::Leapfrog<leapfrog::Py, double, ModelVariant>;

auto get_opts_py(
  const nb::dict parameters,
  const std::vector<int> output_years
) {
  const int t_art_start = nb::cast<int>(parameters["t_ART_start"]);
  const int hts_per_year = nb::cast<int>(parameters["hts_per_year"]);
  const int proj_start_year = nb::cast<int>(parameters["projection_start_year"]);
  const std::string projection_period = nb::cast<std::string>(parameters["projection_period"]);
  return leapfrog::get_opts<double>(
    hts_per_year, t_art_start, projection_period, proj_start_year, output_years
  );
}


template<typename ModelVariant>
auto build_output_py(
  typename LeapfrogPy<ModelVariant>::OutputState output_state,
  const std::vector<int> output_years
) {
  using LF = LeapfrogPy<ModelVariant>;

  nb::dict ret;
  LF::Cfg::build_output(0, output_state, ret, output_years.size());

  return ret;
}

template<typename ModelVariant>
auto build_output_py(
  typename LeapfrogPy<ModelVariant>::State state
) {
  using LF = LeapfrogPy<ModelVariant>;

  nb::dict ret;
  LF::Cfg::build_output_single_year(0, state, ret);

  return ret;
}


template<typename ModelVariant>
nb::dict simulate_model(
  const nb::dict parameters,
  const std::vector<int> output_years
) {
  using LF = LeapfrogPy<ModelVariant>;

  const auto opts = get_opts_py(parameters, output_years);
  const auto pars = LF::Cfg::get_pars(parameters, opts);

  auto state = LF::run_model(pars, opts, output_years);

  return build_output_py<ModelVariant>(state, output_years);
}

template<typename ModelVariant>
nb::dict simulate_model(
  const nb::dict parameters,
  const std::vector<int> output_years,
  nb::dict initial_state_data,
  int simulation_start_year
) {
  using LF = LeapfrogPy<ModelVariant>;

  const auto opts = get_opts_py(parameters, output_years);
  const auto pars = LF::Cfg::get_pars(parameters, opts);

  typename LF::State initial_state = LF::Cfg::get_initial_state(initial_state_data);
  auto state = LF::run_model_from_state(pars, opts, initial_state, simulation_start_year, output_years);

  return build_output_py<ModelVariant>(state, output_years);
}

template<typename ModelVariant>
nb::dict simulate_model(
  const nb::dict parameters,
  nb::dict initial_state_data,
  int simulation_start_year
) {
  using LF = LeapfrogPy<ModelVariant>;

  const auto opts = get_opts_py(parameters, { simulation_start_year + 1 });
  const auto pars = LF::Cfg::get_pars(parameters, opts);

  typename LF::State initial_state = LF::Cfg::get_initial_state(initial_state_data);
  auto state = LF::run_model_single_year(pars, opts, initial_state, simulation_start_year);

  return build_output_py<ModelVariant>(state);
}


std::vector<std::string> list_model_configurations() {
  return std::vector<std::string>{
    "DemographicProjection",
    "HivFullAgeStratification",
    "HivCoarseAgeStratification",
    "ChildModel"
  };
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


nb::dict run_base_model(
  const nb::dict parameters,
  const std::string configuration,
  const std::vector<int> output_years
) {
  return sim_model(configuration, parameters, output_years);
}

nb::dict run_base_model_from_state(
  const nb::dict parameters,
  const std::string configuration,
  const nb::dict initial_state,
  int simulation_start_year,
  const std::vector<int> output_years
) {
  return sim_model(configuration, parameters, output_years, initial_state, simulation_start_year);
}

nb::dict run_base_model_single_year(
  const nb::dict parameters,
  const std::string configuration,
  const nb::dict initial_state,
  int simulation_start_year
) {
  return sim_model(configuration, parameters, initial_state, simulation_start_year);
}

NB_MODULE(_core, m) {
  m.doc() = "Leapfrog python interface";

  m.def("run_base_model", &run_base_model, R"pbdoc(
      Run the leapfrog model.
  )pbdoc");
  m.def("run_base_model_from_state", &run_base_model_from_state, R"pbdoc(
      Run the leapfrog model from an initial state.
  )pbdoc");
  m.def("run_base_model_single_year", &run_base_model_single_year, R"pbdoc(
      Run the leapfrog model from an initial state for a single year.
  )pbdoc");
}
