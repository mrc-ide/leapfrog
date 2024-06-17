#pragma once

#include "generated/parameter_types.hpp"
#include "generated/state_types.hpp"
#include "project_year.hpp"
#include "state_saver.hpp"
#include "frogger.hpp"
#include "utils.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
py::dict build_python_output(const OutputState<ModelVariant, real_type> &state,
                                   const std::vector<int> save_steps) {
  int output_years = save_steps.size();
  constexpr auto ss = StateSpace<ModelVariant>();
  constexpr auto base = ss.base;
  py::array_t<real_type, py::array::f_style> p_total_pop({base.pAG, base.NS, output_years}, state.base.p_total_pop.data());
  py::array_t<real_type, py::array::f_style> births({output_years}, state.base.births.data());
  py::array_t<real_type, py::array::f_style> p_total_pop_natural_deaths({base.pAG, base.NS, output_years}, state.base.p_total_pop_natural_deaths.data());
  py::array_t<real_type, py::array::f_style> p_hiv_pop({base.pAG, base.NS, output_years}, state.base.p_hiv_pop.data());
  py::array_t<real_type, py::array::f_style> p_hiv_pop_natural_deaths({base.pAG, base.NS, output_years}, state.base.p_hiv_pop_natural_deaths.data());
  py::array_t<real_type, py::array::f_style> h_hiv_adult({base.hDS, base.hAG, base.NS, output_years}, state.base.h_hiv_adult.data());
  py::array_t<real_type, py::array::f_style> h_art_adult({base.hTS, base.hDS, base.hAG, base.NS, output_years}, state.base.h_art_adult.data());
  py::array_t<real_type, py::array::f_style> h_hiv_deaths_no_art({base.hDS, base.hAG, base.NS, output_years}, state.base.h_hiv_deaths_no_art.data());
  py::array_t<real_type, py::array::f_style> p_infections({base.pAG, base.NS, output_years}, state.base.p_infections.data());
  py::array_t<real_type, py::array::f_style> h_hiv_deaths_art({base.hTS, base.hDS, base.hAG, base.NS, output_years}, state.base.h_hiv_deaths_art.data());
  py::array_t<real_type, py::array::f_style> h_art_initiation({base.hDS, base.hAG, base.NS, output_years}, state.base.h_art_initiation.data());
  py::array_t<real_type, py::array::f_style> p_hiv_deaths({base.pAG, base.NS, output_years}, state.base.p_hiv_deaths.data());

  py::dict result;
  result["p_total_pop"] = p_total_pop;
  result["births"] = births;
  result["p_total_pop_natural_deaths"] = p_total_pop_natural_deaths;
  result["p_hiv_pop"] = p_hiv_pop;
  result["p_hiv_pop_natural_deaths"] = p_hiv_pop_natural_deaths;
  result["h_hiv_adult"] = h_hiv_adult;
  result["h_art_adult"] = h_art_adult;
  result["h_hiv_deaths_no_art"] = h_hiv_deaths_no_art;
  result["p_infections"] = p_infections;
  result["h_hiv_deaths_art"] = h_hiv_deaths_art;
  result["h_art_initiation"] = h_art_initiation;
  result["p_hiv_deaths"] = p_hiv_deaths;

  if constexpr (ModelVariant::run_child_model) {
    constexpr auto children = ss.children;
    constexpr auto base = ss.base;
    py::array_t<real_type, py::array::f_style> hc1_hiv_pop({children.hc1DS, children.hcTT, children.hc1AG, base.NS, output_years}, state.children.hc1_hiv_pop.data());
    py::array_t<real_type, py::array::f_style> hc2_hiv_pop({children.hc2DS, children.hcTT, children.hc2AG, base.NS, output_years}, state.children.hc2_hiv_pop.data());
    py::array_t<real_type, py::array::f_style> hc1_art_pop({base.hTS, children.hc1DS, children.hc1AG, base.NS, output_years}, state.children.hc1_art_pop.data());
    py::array_t<real_type, py::array::f_style> hc2_art_pop({base.hTS, children.hc2DS, children.hc2AG, base.NS, output_years}, state.children.hc2_art_pop.data());
    py::array_t<real_type, py::array::f_style> hc1_noart_aids_deaths({children.hc1DS, children.hcTT, children.hc1AG, base.NS, output_years}, state.children.hc1_noart_aids_deaths.data());
    py::array_t<real_type, py::array::f_style> hc2_noart_aids_deaths({children.hc2DS, children.hcTT, children.hc2AG, base.NS, output_years}, state.children.hc2_noart_aids_deaths.data());
    py::array_t<real_type, py::array::f_style> hc1_art_aids_deaths({base.hTS, children.hc1DS, children.hc1AG, base.NS, output_years}, state.children.hc1_art_aids_deaths.data());
    py::array_t<real_type, py::array::f_style> hc2_art_aids_deaths({base.hTS, children.hc2DS, children.hc2AG, base.NS, output_years}, state.children.hc2_art_aids_deaths.data());
    py::array_t<real_type, py::array::f_style> hiv_births({output_years}, state.children.hiv_births.data());
    py::array_t<real_type, py::array::f_style> hc_art_total({4, output_years}, state.children.hc_art_total.data());
    py::array_t<real_type, py::array::f_style> hc_art_init({4, output_years}, state.children.hc_art_init.data());
    py::array_t<real_type, py::array::f_style> hc_art_need_init({children.hc1DS, children.hcTT, 15, base.NS, output_years}, state.children.hc_art_need_init.data());
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
    result["hc_art_total"] = hc_art_total;
    result["hc_art_init"] = hc_art_init;
    result["hc_art_need_init"] = hc_art_need_init;
    result["ctx_need"] = ctx_need;
    result["ctx_mean"] = ctx_mean;
  }

  return result;
}

template<typename ModelVariant, typename real_type>
leapfrog::Parameters<ModelVariant, real_type> setup_model_params_py(const py::dict &params,
                                                                    const leapfrog::Options<real_type> &options) {
  constexpr auto ss = leapfrog::StateSpace<ModelVariant>();

  constexpr auto base = ss.base;
  const leapfrog::TensorMap2<real_type> base_pop = params["basepop"].cast<leapfrog::TensorMap2<real_type>>();
  const leapfrog::TensorMap3<real_type> survival_probability = params["Sx"].cast<leapfrog::TensorMap3<real_type>>();
  const leapfrog::TensorMap3<real_type> net_migration = params["netmigr_adj"].cast<leapfrog::TensorMap3<real_type>>();
  const leapfrog::TensorMap2<real_type> age_specific_fertility_rate = params["asfr"].cast<leapfrog::TensorMap2<real_type>>();
  const leapfrog::TensorMap2<real_type> births_sex_prop = params["births_sex_prop"].cast<leapfrog::TensorMap2<real_type>>();
  const leapfrog::TensorMap1<real_type> total_rate = params["incidinput"].cast<leapfrog::TensorMap1<real_type>>();
  const leapfrog::TensorMap3<real_type> relative_risk_age = params["incrr_age"].cast<leapfrog::TensorMap3<real_type>>();
  const leapfrog::TensorMap1<real_type> relative_risk_sex = params["incrr_sex"].cast<leapfrog::TensorMap1<real_type>>();
  const leapfrog::TensorMap3<real_type> cd4_mortality = params["cd4_mort"].cast<leapfrog::TensorMap3<real_type>>();
  const leapfrog::TensorMap3<real_type> cd4_progression = params["cd4_prog"].cast<leapfrog::TensorMap3<real_type>>();
  const leapfrog::Tensor1<int> idx_hm_elig = convert_0_based<1>(params["artcd4elig_idx"].cast<leapfrog::TensorMap1<int>>());
  const leapfrog::TensorMap3<real_type> cd4_initial_distribution = params["cd4_initdist"].cast<leapfrog::TensorMap3<real_type>>();
  const leapfrog::TensorMap4<real_type> mortality = params["art_mort"].cast<leapfrog::TensorMap4<real_type>>();
  const leapfrog::TensorMap2<real_type> mortaility_time_rate_ratio = params["artmx_timerr"].cast<leapfrog::TensorMap2<real_type>>();
  const leapfrog::TensorMap1<real_type> dropout = params["art_dropout"].cast<leapfrog::TensorMap1<real_type>>();
  const leapfrog::TensorMap2<real_type> adults_on_art = params["art15plus_num"].cast<leapfrog::TensorMap2<real_type>>();
  const leapfrog::TensorMap2<int> adults_on_art_is_percent = params["art15plus_isperc"].cast<leapfrog::TensorMap2<int>>();
  const int scale_cd4_mortality = params["scale_cd4_mort"].cast<int>();
  const real_type initiation_mortality_weight = params["art_alloc_mxweight"].cast<real_type>();
  leapfrog::Tensor1<real_type> h_art_stage_dur(base.hTS - 1);
  h_art_stage_dur.setConstant(0.5);

  const leapfrog::Demography<real_type> demography_params = {
    base_pop,
    survival_probability,
    net_migration,
    age_specific_fertility_rate,
    births_sex_prop,
  };
  const leapfrog::Incidence<real_type> incidence_params = {
    total_rate,
    relative_risk_age,
    relative_risk_sex,
  };
  const leapfrog::NaturalHistory<real_type> natural_history_params = {
    cd4_mortality,
    cd4_progression,
    cd4_initial_distribution,
    scale_cd4_mortality,
  };
  const leapfrog::Art<real_type> art_params = {
    idx_hm_elig,
    mortality,
    mortaility_time_rate_ratio,
    dropout,
    adults_on_art,
    adults_on_art_is_percent,
    h_art_stage_dur,
    initiation_mortality_weight,
  };


  const leapfrog::BaseModelParameters<real_type> base_model_params = {
    options,
    demography_params,
    incidence_params,
    natural_history_params,
    art_params
  };

  if constexpr (ModelVariant::run_child_model) {
    constexpr auto children = ss.children;
    constexpr auto base = ss.base;
    const leapfrog::TensorMap1<real_type> hc_nosocomial = params["paed_incid_input"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::TensorMap1<real_type> hc1_cd4_dist = params["paed_cd4_dist"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::TensorMap2<real_type> hc_cd4_transition = params["paed_cd4_transition"].cast<leapfrog::TensorMap2<real_type>>();
    const leapfrog::TensorMap3<real_type> hc1_cd4_mort = params["paed_cd4_mort"].cast<leapfrog::TensorMap3<real_type>>();
    const leapfrog::TensorMap3<real_type> hc2_cd4_mort = params["adol_cd4_mort"].cast<leapfrog::TensorMap3<real_type>>();
    const leapfrog::TensorMap1<real_type> hc1_cd4_prog = params["paed_cd4_prog"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::TensorMap1<real_type> hc2_cd4_prog = params["adol_cd4_prog"].cast<leapfrog::TensorMap1<real_type>>();
    const real_type ctx_effect = params["ctx_effect"].cast<real_type>();
    const leapfrog::TensorMap1<real_type> ctx_val = params["ctx_val"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::TensorMap1<int> hc_art_elig_age = params["paed_art_elig_age"].cast<leapfrog::TensorMap1<int>>();
    const leapfrog::Tensor2<real_type> hc_art_elig_cd4 = convert_0_based<2>(params["paed_art_elig_cd4"].cast<leapfrog::TensorMap2<real_type>>());
    const leapfrog::TensorMap3<real_type> hc_art_mort_rr = params["mort_art_rr"].cast<leapfrog::TensorMap3<real_type>>();
    const leapfrog::TensorMap3<real_type> hc1_art_mort = params["paed_art_mort"].cast<leapfrog::TensorMap3<real_type>>();
    const leapfrog::TensorMap3<real_type> hc2_art_mort = params["adol_art_mort"].cast<leapfrog::TensorMap3<real_type>>();
    const leapfrog::TensorMap1<int> hc_art_isperc = params["artpaeds_isperc"].cast<leapfrog::TensorMap1<int>>();
    const leapfrog::TensorMap2<real_type> hc_art_val = params["paed_art_val"].cast<leapfrog::TensorMap2<real_type>>();
    const leapfrog::TensorMap2<real_type> hc_art_init_dist = params["init_art_dist"].cast<leapfrog::TensorMap2<real_type>>();
    const leapfrog::TensorMap2<real_type> adult_cd4_dist = params["adult_cd4_dist"].cast<leapfrog::TensorMap2<real_type>>();
    const leapfrog::TensorMap1<real_type> fert_mult_by_age = params["fert_mult_by_age"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::TensorMap1<real_type> fert_mult_off_art = params["fert_mult_offart"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::TensorMap1<real_type> fert_mult_on_art = params["fert_mult_onart"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::TensorMap1<real_type> total_fertility_rate = params["tfr"].cast<leapfrog::TensorMap1<real_type>>();
    const real_type local_adj_factor = params["laf"].cast<real_type>();
    const leapfrog::TensorMap2<real_type> PMTCT = params["pmtct"].cast<leapfrog::TensorMap2<real_type>>();
    const leapfrog::TensorMap2<real_type> vertical_transmission_rate = params["mtct"].cast<leapfrog::TensorMap2<real_type>>();
    const leapfrog::TensorMap3<real_type> PMTCT_transmission_rate = params["pmtct_mtct"].cast<leapfrog::TensorMap3<real_type>>();
    const leapfrog::TensorMap2<real_type> PMTCT_dropout = params["pmtct_dropout"].cast<leapfrog::TensorMap2<real_type>>();
    const leapfrog::TensorMap1<int> PMTCT_input_is_percent = params["pmtct_input_isperc"].cast<leapfrog::TensorMap1<int>>();
    const leapfrog::TensorMap2<real_type> breastfeeding_duration_art = params["bf_duration_art"].cast<leapfrog::TensorMap2<real_type>>();
    const leapfrog::TensorMap2<real_type> breastfeeding_duration_no_art = params["bf_duration_no_art"].cast<leapfrog::TensorMap2<real_type>>();
    const leapfrog::TensorMap1<real_type> mat_hiv_births = params["mat_hiv_births"].cast<leapfrog::TensorMap2<real_type>>();
    const leapfrog::TensorMap1<int> mat_prev_input = params["mat_prev_input"].cast<leapfrog::TensorMap1<int>>();
    const leapfrog::TensorMap1<real_type> prop_lt200 = params["prop_lt200"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::TensorMap1<real_type> prop_gte350 = params["prop_gte350"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::TensorMap1<real_type> incrate = params["incrate"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::TensorMap1<int> ctx_val_is_percent = params["ctx_val_ispercent"].cast<leapfrog::TensorMap1<int>>();
    const leapfrog::TensorMap1<int> hc_art_is_age_spec = params["paed_art_age_spec"].cast<leapfrog::TensorMap1<int>>();
    const leapfrog::TensorMap1<real_type> hc_age_coarse = params["hc_age_coarse"].cast<leapfrog::TensorMap1<real_type>>();
    const leapfrog::Children<real_type> children_params = {
      hc_nosocomial,
      hc1_cd4_dist,
      hc_cd4_transition,
      hc1_cd4_mort,
      hc2_cd4_mort,
      hc1_cd4_prog,
      hc2_cd4_prog,
      ctx_effect,
      ctx_val,
      hc_art_elig_age,
      hc_art_elig_cd4,
      hc_art_mort_rr,
      hc1_art_mort,
      hc2_art_mort,
      hc_art_isperc,
      hc_art_val,
      hc_art_init_dist,
      adult_cd4_dist,
      fert_mult_by_age,
      fert_mult_off_art,
      fert_mult_on_art,
      total_fertility_rate,
      local_adj_factor,
      PMTCT,
      vertical_transmission_rate,
      PMTCT_transmission_rate,
      PMTCT_dropout,
      PMTCT_input_is_percent,
      breastfeeding_duration_art,
      breastfeeding_duration_no_art,
      mat_hiv_births,
      mat_prev_input,
      prop_lt200,
      prop_gte350,
      incrate,
      ctx_val_is_percent,
      hc_art_is_age_spec,
      hc_age_coarse,
    };
    const leapfrog::ChildModelParameters<ModelVariant, real_type> child_model_params = {
      children_params
    };
    return leapfrog::Parameters<ModelVariant, real_type> {
      base_model_params,
      child_model_params
    };
  } else {
    return leapfrog::Parameters<ModelVariant, real_type> {
      base_model_params
    };
  }
}

} // namespace internal

template<typename ModelVariant, typename real_type, bool OwnedData>
void set_initial(State<ModelVariant, real_type, OwnedData> &state,
                 const py::dict &parameters) {

  constexpr auto ss = leapfrog::StateSpace<ModelVariant>();
  const leapfrog::Options<double> opts = {
    0, // HIV steps not needed here, as not running the model
    parameters["t_ART_start"].cast<int>() - 1,
    ss.base.hAG
  };

  const auto params = internal::setup_model_params_py<ModelVariant, double>(parameters, opts);

  set_initial_state(state, params);
}

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
template<typename ModelVariant, typename real_type, bool OwnedData>
void project_single_year(int time_step,
                         const int hiv_steps,
                         const py::dict &parameters,
                         const State<ModelVariant, real_type, OwnedData> &state_curr,
                         State<ModelVariant, real_type, OwnedData> &state_next) {

  constexpr auto ss = leapfrog::StateSpace<ModelVariant>();
  const leapfrog::Options<double> opts = {
    hiv_steps,
    parameters["t_ART_start"].cast<int>() - 1,
    ss.base.hAG
  };

  const auto params = internal::setup_model_params_py<ModelVariant, real_type>(parameters, opts);

  internal::IntermediateData<ModelVariant, real_type> intermediate(params.base.options.hAG_15plus);
  intermediate.reset();

  internal::project_year(time_step, params, state_curr, state_next, intermediate);
}


template<typename ModelVariant, typename real_type>
py::dict simulate_model(const py::dict &parameters,
                        const int proj_years,
                        const int hiv_steps,
                        const std::vector<int> save_steps) {

  constexpr auto ss = leapfrog::StateSpace<ModelVariant>();
  const leapfrog::Options<double> opts = {
    hiv_steps,
    parameters["t_ART_start"].cast<int>() - 1,
    ss.base.hAG
  };

  const auto params = internal::setup_model_params_py<ModelVariant, double>(parameters, opts);

  auto state = run_model<ModelVariant, real_type>(proj_years, save_steps, params);

  auto ret = internal::build_python_output<ModelVariant, real_type>(state, save_steps);

  return ret;
}

} // namespace leapfrog
