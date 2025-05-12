#include "../inst/include/frogger.hpp"
#include "../inst/include/generated/tensor_types.hpp"
#include "../inst/include/generated/c_interface/c_types.hpp"
#include "../inst/include/generated/c_interface/c_adapters.hpp"

#include <windows.h>
#include <numeric>
#include <stdexcept>

#define DllExport extern "C"
#define EXPORT comment(linker, "/EXPORT:" __FUNCTION__ "=" __FUNCDNAME__)

typedef void (WINAPI *CallbackFunction)(const char*);

template <typename ModelVariant>
HRESULT fit_model(leapfrog::internal::COptions &options,
                  leapfrog::internal::CParams<double> &data,
                  leapfrog::internal::COutput<double> &out,
                  CallbackFunction error_handler) {

  using LF = leapfrog::Leapfrog<leapfrog::C, double, ModelVariant>;

  try {
    std::vector<int> output_years((options.proj_end_year - options.proj_start_year) + 1);
    std::iota(output_years.begin(), output_years.end(), options.proj_start_year);

    const leapfrog::Options<double> opts =  {
      10,
      34,
      leapfrog::internal::BaseSS::PROJPERIOD_CALENDAR,
      options.proj_start_year,
      options.proj_end_year
    };
    const auto pars = LF::Cfg::get_pars(opts, data);

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

DllExport HRESULT WINAPI run_dp(leapfrog::internal::COptions &options,
                                leapfrog::internal::CParams<double> &data,
                                leapfrog::internal::COutput<double> &out,
                                CallbackFunction error_handler) {
  #pragma EXPORT
  return fit_model<leapfrog::DemographicProjection>(options, data, out, error_handler);
}

DllExport HRESULT WINAPI run_aim(leapfrog::internal::COptions &options,
                                 leapfrog::internal::CParams<double> &data,
                                 leapfrog::internal::COutput<double> &out,
                                 CallbackFunction error_handler) {
  #pragma EXPORT
  return fit_model<leapfrog::ChildModel>(options, data, out, error_handler);
}

