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

DllExport HRESULT WINAPI run_dp(leapfrog::internal::DpParams& data, leapfrog::internal::DpOut& out, CallbackFunction logger) {
  #pragma EXPORT

  using LF = leapfrog::Leapfrog<leapfrog::C, double, leapfrog::DemographicProjection>;

  try {
    std::vector<int> output_years(61);
    std::iota(output_years.begin(), output_years.end(), 1970);

    const auto opts = leapfrog::get_opts<double>(10, 34, std::string_view{"calendar"}, 1970, output_years);
    const auto pars = LF::Cfg::get_pars(data, opts);

    auto state = LF::run_model(pars, opts, output_years);
    LF::Cfg::build_output(0, state, out);
  } catch (const std::invalid_argument& e) {
    logger(e.what());
    return E_INVALIDARG;
  } catch (...) {
    logger("Caught unhandled exception");
    return E_FAIL;
  }

  return S_OK;
}
