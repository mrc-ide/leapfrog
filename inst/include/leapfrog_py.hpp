#pragma once

#include "generated/parameter_types.hpp"
#include "generated/state_types.hpp"
#include "project_year.hpp"
#include "state_saver.hpp"
#include "frogger.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace leapfrog {

/**
 * @brief Run a simulation model over a specified number of time steps.
 *
 * @tparam ModelVariant The variant of the model to be run, used for compile time switching.
 * @tparam real_type The data type used for real numbers in the simulation, usually a double.
 * @param time_step The time step to run, used in index parameters
 * @param pars The parameters required for running the simulation, read only.
 * @param state_curr The current state of the model, read only.
 * @param state_next The next state of the model.
 * @return None, updates state_next in place
 */
template<typename ModelVariant, typename real_type>
void project_single_year(int time_step,
                         const Parameters<ModelVariant, real_type> &pars,
                         const State<ModelVariant, real_type> &state_curr,
                         State<ModelVariant, real_type> &state_next) {

  internal::IntermediateData<ModelVariant, real_type> intermediate(pars.options.hAG_15plus);
  intermediate.reset();

  internal::project_year(time_step, pars, state_curr, state_next, intermediate);
}

namespace internal {

template<typename ModelVariant, typename real_type>
py::dict build_python_output(const OutputState<ModelVariant, real_type> &state,
                             const std::vector<int> save_steps) {
  int output_years = save_steps.size();
  constexpr auto ss = StateSpace<ModelVariant>();
  constexpr auto dp = ss.dp;
  py::array_t<real_type, py::array::f_style> p_total_pop({dp.pAG, dp.NS, output_years}, state.dp.p_total_pop.data());
  py::array_t<real_type, py::array::f_style> births({output_years}, state.dp.births.data());
  py::array_t<real_type, py::array::f_style> p_total_pop_natural_deaths({dp.pAG, dp.NS, output_years}, state.dp.p_total_pop_natural_deaths.data());

  py::dict result;
  result["p_total_pop"] = p_total_pop;
  result["births"] = births;
  result["p_total_pop_natural_deaths"] = p_total_pop_natural_deaths;

  if constexpr (ModelVariant::run_hiv_simulation) {
    constexpr auto dp = ss.dp;
    constexpr auto hiv = ss.hiv;
    py::array_t<real_type, py::array::f_style> p_hiv_pop({dp.pAG, dp.NS, output_years}, state.hiv.p_hiv_pop.data());
    py::array_t<real_type, py::array::f_style> p_hiv_pop_natural_deaths({dp.pAG, dp.NS, output_years}, state.hiv.p_hiv_pop_natural_deaths.data());
    py::array_t<real_type, py::array::f_style> h_hiv_adult({hiv.hDS, hiv.hAG, dp.NS, output_years}, state.hiv.h_hiv_adult.data());
    py::array_t<real_type, py::array::f_style> h_art_adult({hiv.hTS, hiv.hDS, hiv.hAG, dp.NS, output_years}, state.hiv.h_art_adult.data());
    py::array_t<real_type, py::array::f_style> h_hiv_deaths_no_art({hiv.hDS, hiv.hAG, dp.NS, output_years}, state.hiv.h_hiv_deaths_no_art.data());
    py::array_t<real_type, py::array::f_style> p_infections({dp.pAG, dp.NS, output_years}, state.hiv.p_infections.data());
    py::array_t<real_type, py::array::f_style> h_hiv_deaths_art({hiv.hTS, hiv.hDS, hiv.hAG, dp.NS, output_years}, state.hiv.h_hiv_deaths_art.data());
    py::array_t<real_type, py::array::f_style> h_art_initiation({hiv.hDS, hiv.hAG, dp.NS, output_years}, state.hiv.h_art_initiation.data());
    py::array_t<real_type, py::array::f_style> p_hiv_deaths({dp.pAG, dp.NS, output_years}, state.hiv.p_hiv_deaths.data());

    result["p_hiv_pop"] = p_hiv_pop;
    result["p_hiv_pop_natural_deaths"] = p_hiv_pop_natural_deaths;
    result["h_hiv_adult"] = h_hiv_adult;
    result["h_art_adult"] = h_art_adult;
    result["h_hiv_deaths_no_art"] = h_hiv_deaths_no_art;
    result["p_infections"] = p_infections;
    result["h_hiv_deaths_art"] = h_hiv_deaths_art;
    result["h_art_initiation"] = h_art_initiation;
    result["p_hiv_deaths"] = p_hiv_deaths;
  }

  if constexpr (ModelVariant::run_child_model) {
    constexpr auto children = ss.children;
    constexpr auto dp = ss.dp;
    constexpr auto hiv = ss.hiv;

    py::array_t<real_type, py::array::f_style> hc1_hiv_pop({children.hc1DS, children.hcTT, children.hc1AG, dp.NS, output_years}, state.children.hc1_hiv_pop.data());
    py::array_t<real_type, py::array::f_style> hc2_hiv_pop({children.hc2DS, children.hcTT, children.hc2AG, dp.NS, output_years}, state.children.hc2_hiv_pop.data());
    py::array_t<real_type, py::array::f_style> hc1_art_pop({hiv.hTS, children.hc1DS, children.hc1AG, dp.NS, output_years}, state.children.hc1_art_pop.data());
    py::array_t<real_type, py::array::f_style> hc2_art_pop({hiv.hTS, children.hc2DS, children.hc2AG, dp.NS, output_years}, state.children.hc2_art_pop.data());
    py::array_t<real_type, py::array::f_style> hc1_noart_aids_deaths({children.hc1DS, children.hcTT, children.hc1AG, dp.NS, output_years}, state.children.hc1_noart_aids_deaths.data());
    py::array_t<real_type, py::array::f_style> hc2_noart_aids_deaths({children.hc2DS, children.hcTT, children.hc2AG, dp.NS, output_years}, state.children.hc2_noart_aids_deaths.data());
    py::array_t<real_type, py::array::f_style> hc1_art_aids_deaths({hiv.hTS, children.hc1DS, children.hc1AG, dp.NS, output_years}, state.children.hc1_art_aids_deaths.data());
    py::array_t<real_type, py::array::f_style> hc2_art_aids_deaths({hiv.hTS, children.hc2DS, children.hc2AG, dp.NS, output_years}, state.children.hc2_art_aids_deaths.data());
    py::array_t<real_type, py::array::f_style> hiv_births({output_years}, state.children.hiv_births.data());
    py::array_t<real_type, py::array::f_style> hc_art_init({children.hcAG_coarse, output_years}, state.children.hc_art_init.data());
    py::array_t<real_type, py::array::f_style> hc_art_need_init({children.hc1DS, children.hcTT, children.hcAG_end, dp.NS, output_years}, state.children.hc_art_need_init.data());
    py::array_t<real_type, py::array::f_style> ctx_need({output_years}, state.children.ctx_need.data());
    py::array_t<real_type, py::array::f_style> ctx_mean({output_years}, state.children.ctx_mean.data());

    result["hc1_hiv_pop"] = hc1_hiv_pop;
    result["hc2_hiv_pop"] = hc2_hiv_pop;
    result["hc1_art_pop"] = hc1_art_pop;
    result["hc2_art_pop"] = hc2_art_pop;
    result["hc1_noart_aids_deaths"] = hc1_noart_aids_deaths;
    result["hc2_noart_aids_deaths"] = hc2_noart_aids_deaths;
    result["hc1_art_aids_deaths"] = hc1_art_aids_deaths;
    result["hc2_art_aids_deaths"] = hc2_art_aids_deaths;
    result["hiv_births"] = hiv_births;
    result["hc_art_init"] = hc_art_init;
    result["hc_art_need_init"] = hc_art_need_init;
    result["ctx_need"] = ctx_need;
    result["ctx_mean"] = ctx_mean;
  }

  return result;
}

} // namespace internal

template<typename ModelVariant, typename real_type>
py::dict simulate_model(const leapfrog::Parameters<ModelVariant, real_type> &params,
                        const int proj_years,
                        const std::vector<int> save_steps) {

  auto state = run_model<ModelVariant, real_type>(proj_years, save_steps, params);

  auto ret = internal::build_python_output<ModelVariant, real_type>(state, save_steps);

  return ret;
}

} // namespace leapfrog
