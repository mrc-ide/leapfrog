#include "../r-package/inst/include/frogger.hpp"
#include "../r-package/inst/include/generated/tensor_types.hpp"
#include "../r-package/inst/include/generated/c_interface/c_types.hpp"
#include "../r-package/inst/include/generated/c_interface/c_adapters.hpp"

#include <windows.h>
#include <numeric>
#include <stdexcept>

#define DllExport extern "C"
#define EXPORT comment(linker, "/EXPORT:" __FUNCTION__ "=" __FUNCDNAME__)

typedef void (WINAPI *CallbackFunction)(const char*);

template <typename ModelVariant>
HRESULT fit_model(leapfrog::internal::CParams<double> &data,
                  leapfrog::internal::COptions &options,
                  leapfrog::internal::CState<double> &out,
                  CallbackFunction error_handler) {

  using LF = leapfrog::Leapfrog<leapfrog::C, double, ModelVariant>;

  try {
    std::vector<int> output_years((options.proj_end_year - options.proj_start_year) + 1);
    std::iota(output_years.begin(), output_years.end(), options.proj_start_year);

    const leapfrog::Options<double> opts =  {
      10,
      options.ts_art_start,
      options.proj_period,
      options.proj_start_year,
      options.proj_end_year
    };
    const auto pars = LF::Cfg::get_pars(data, opts);

    auto state = LF::run_model(pars, opts, output_years);
    LF::Cfg::build_output(0, state, out);
  } catch (const std::invalid_argument& e) {
    error_handler(e.what());
    return E_INVALIDARG;
  } catch (...) {
    error_handler("Caught unhandled exception");
    return E_FAIL;
  }

  return S_OK;
}

DllExport HRESULT WINAPI run_dp(leapfrog::internal::CParams<double> &data,
  leapfrog::internal::COptions &options,
  leapfrog::internal::CState<double> &out,
  CallbackFunction error_handler) {
#pragma EXPORT
return fit_model<leapfrog::DemographicProjection>(data, options, out, error_handler);
}

DllExport HRESULT WINAPI run_aim(leapfrog::internal::CParams<double> &data,
   leapfrog::internal::COptions &options,
   leapfrog::internal::CState<double> &out,
   CallbackFunction error_handler) {
#pragma EXPORT
return fit_model<leapfrog::ChildModel>(data, options, out, error_handler);
}

template <typename ModelVariant>
HRESULT fit_model_single_year(leapfrog::internal::CParams<double> &data,
                              leapfrog::internal::COptions &options,
                              leapfrog::internal::CState<double> &initial_state_data,
                              int simulation_start_year,
                              leapfrog::internal::CState<double> &out,
                              CallbackFunction error_handler) {

  using LF = leapfrog::Leapfrog<leapfrog::C, double, ModelVariant>;

  try {
    const leapfrog::Options<double> opts =  {
      10,
      options.ts_art_start,
      options.proj_period,
      options.proj_start_year,
      options.proj_end_year
    };
    const auto pars = LF::Cfg::get_pars(data, opts);
    const auto initial_state = LF::Cfg::get_initial_state(initial_state_data);

    auto state = LF::run_model_single_year(pars, opts, initial_state, simulation_start_year);
    LF::Cfg::build_output_single_year(0, state, out);
  } catch (const std::invalid_argument& e) {
    error_handler(e.what());
    return E_INVALIDARG;
  } catch (...) {
    error_handler("Caught unhandled exception");
    return E_FAIL;
  }

  return S_OK;
}

DllExport HRESULT WINAPI run_dp_single_year(leapfrog::internal::CParams<double> &data,
                                            leapfrog::internal::COptions &options,
                                            leapfrog::internal::CState<double> &initial_state,
                                            int simulation_start_year,
                                            leapfrog::internal::CState<double> &out,
                                            CallbackFunction error_handler) {
  #pragma EXPORT
  return fit_model_single_year<leapfrog::DemographicProjection>(data, options, initial_state, simulation_start_year, out, error_handler);
}

DllExport HRESULT WINAPI run_aim_single_year(leapfrog::internal::CParams<double> &data,
                                             leapfrog::internal::COptions &options,
                                             leapfrog::internal::CState<double> &initial_state,
                                             int simulation_start_year,
                                             leapfrog::internal::CState<double> &out,
                                             CallbackFunction error_handler) {
  #pragma EXPORT
  return fit_model_single_year<leapfrog::ChildModel>(data, options, initial_state, simulation_start_year, out, error_handler);
}

