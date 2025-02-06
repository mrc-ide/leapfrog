#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen/tensor.h>
#include "../../inst/include/generated/parameter_types.hpp"
#include "../../inst/include/generated/state_types.hpp"
#include "../../inst/include/leapfrog_py.hpp"
#include "../../inst/include/frogger.hpp"

namespace py = pybind11;

namespace {
using Eigen::Sizes;
}

PYBIND11_MODULE(leapfrog, m) {
    m.doc() = "pybind11 leapfrog plugin";

    py::class_<leapfrog::Demography<double>>(m, "Demography")
        .def(py::init<const leapfrog::TensorMap2<double>&,
                      const leapfrog::TensorMap3<double>&,
                      const leapfrog::TensorMap3<double>&,
                      const leapfrog::TensorMap2<double>&,
                      const leapfrog::TensorMap2<double>&>())
        .def_readonly("base_pop", &leapfrog::Demography<double>::base_pop)
        .def_readonly("survival_probability", &leapfrog::Demography<double>::survival_probability)
        .def_readonly("net_migration", &leapfrog::Demography<double>::net_migration)
        .def_readonly("age_specific_fertility_rate", &leapfrog::Demography<double>::age_specific_fertility_rate)
        .def_readonly("births_sex_prop", &leapfrog::Demography<double>::births_sex_prop);

    py::class_<leapfrog::Incidence<double>>(m, "Incidence")
        .def(py::init<const leapfrog::TensorMap1<double>&,
                      const leapfrog::TensorMap3<double>&,
                      const leapfrog::TensorMap1<double>&>())
        .def_readonly("total_rate", &leapfrog::Incidence<double>::total_rate)
        .def_readonly("relative_risk_age", &leapfrog::Incidence<double>::relative_risk_age)
        .def_readonly("relative_risk_sex", &leapfrog::Incidence<double>::relative_risk_sex);

    py::class_<leapfrog::NaturalHistory<double>>(m, "NaturalHistory")
        .def(py::init<const leapfrog::TensorMap3<double>&,
                      const leapfrog::TensorMap3<double>&,
                      const leapfrog::TensorMap3<double>&,
                      int>())
        .def_readonly("cd4_mortality", &leapfrog::NaturalHistory<double>::cd4_mortality)
        .def_readonly("cd4_progression", &leapfrog::NaturalHistory<double>::cd4_progression)
        .def_readonly("cd4_initial_distribution", &leapfrog::NaturalHistory<double>::cd4_initial_distribution)
        .def_readonly("scale_cd4_mortality", &leapfrog::NaturalHistory<double>::scale_cd4_mortality);

    py::class_<leapfrog::Art<double>>(m, "Art")
        .def(py::init<const leapfrog::Tensor1<int>&,
                      const leapfrog::TensorMap4<double>&,
                      const leapfrog::TensorMap2<double>&,
                      int,
                      const leapfrog::TensorMap1<double>&,
                      const leapfrog::TensorMap2<double>&,
                      const leapfrog::TensorMap2<int>&,
                      const leapfrog::Tensor1<double>&,
                      double>())
        .def_readonly("idx_hm_elig", &leapfrog::Art<double>::idx_hm_elig)
        .def_readonly("mortality", &leapfrog::Art<double>::mortality)
        .def_readonly("mortality_time_rate_ratio", &leapfrog::Art<double>::mortality_time_rate_ratio)
        .def_readonly("dropout_recover_cd4", &leapfrog::Art<double>::dropout_recover_cd4)
        .def_readonly("dropout_rate", &leapfrog::Art<double>::dropout_rate)
        .def_readonly("adults_on_art", &leapfrog::Art<double>::adults_on_art)
        .def_readonly("adults_on_art_is_percent", &leapfrog::Art<double>::adults_on_art_is_percent)
        .def_readonly("h_art_stage_dur", &leapfrog::Art<double>::h_art_stage_dur)
        .def_readonly("initiation_mortality_weight", &leapfrog::Art<double>::initiation_mortality_weight);

    py::class_<leapfrog::Children<double>>(m, "Children")
            .def(py::init<const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap3<double>&,
                          const leapfrog::TensorMap3<double>&,
                          const leapfrog::TensorMap3<double>&,
                          const leapfrog::TensorMap3<double>&,
                          const double&,
                          const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap1<int>&,
                          const leapfrog::Tensor2<double>&,
                          const leapfrog::TensorMap3<double>&,
                          const leapfrog::TensorMap3<double>&,
                          const leapfrog::TensorMap3<double>&,
                          const leapfrog::TensorMap1<int>&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap1<double>&,
                          const double&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap3<double>&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap1<int>&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap1<int>&,

                          // Nans as first 6 numbers here
                          const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap1<double>&,

                          const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap1<int>&,
                          const leapfrog::TensorMap1<int>&,
                          const leapfrog::TensorMap1<double>&,

                          // Should be 3D
                          const leapfrog::TensorMap2<double>&,

                          const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap1<double>&,
                          const leapfrog::TensorMap1<int>&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap2<double>&,
                          const leapfrog::TensorMap1<double>&>())
        .def_readonly("hc_nosocomial", &leapfrog::Children<double>::hc_nosocomial)
        .def_readonly("hc1_cd4_dist", &leapfrog::Children<double>::hc1_cd4_dist)
        .def_readonly("hc_cd4_transition", &leapfrog::Children<double>::hc_cd4_transition)
        .def_readonly("hc1_cd4_mort", &leapfrog::Children<double>::hc1_cd4_mort)
        .def_readonly("hc2_cd4_mort", &leapfrog::Children<double>::hc2_cd4_mort)
        .def_readonly("hc1_cd4_prog", &leapfrog::Children<double>::hc1_cd4_prog)
        .def_readonly("hc2_cd4_prog", &leapfrog::Children<double>::hc2_cd4_prog)
        .def_readonly("ctx_effect", &leapfrog::Children<double>::ctx_effect)
        .def_readonly("ctx_val", &leapfrog::Children<double>::ctx_val)
        .def_readonly("hc_art_elig_age", &leapfrog::Children<double>::hc_art_elig_age)
        .def_readonly("hc_art_elig_cd4", &leapfrog::Children<double>::hc_art_elig_cd4)
        .def_readonly("hc_art_mort_rr", &leapfrog::Children<double>::hc_art_mort_rr)
        .def_readonly("hc1_art_mort", &leapfrog::Children<double>::hc1_art_mort)
        .def_readonly("hc2_art_mort", &leapfrog::Children<double>::hc2_art_mort)
        .def_readonly("hc_art_isperc", &leapfrog::Children<double>::hc_art_isperc)
        .def_readonly("hc_art_val", &leapfrog::Children<double>::hc_art_val)
        .def_readonly("hc_art_init_dist", &leapfrog::Children<double>::hc_art_init_dist)
        .def_readonly("adult_cd4_dist", &leapfrog::Children<double>::adult_cd4_dist)
        .def_readonly("fert_mult_by_age", &leapfrog::Children<double>::fert_mult_by_age)
        .def_readonly("fert_mult_off_art", &leapfrog::Children<double>::fert_mult_off_art)
        .def_readonly("fert_mult_on_art", &leapfrog::Children<double>::fert_mult_on_art)
        .def_readonly("total_fertility_rate", &leapfrog::Children<double>::total_fertility_rate)
        .def_readonly("local_adj_factor", &leapfrog::Children<double>::local_adj_factor)
        .def_readonly("PMTCT", &leapfrog::Children<double>::PMTCT)
        .def_readonly("vertical_transmission_rate", &leapfrog::Children<double>::vertical_transmission_rate)
        .def_readonly("PMTCT_transmission_rate", &leapfrog::Children<double>::PMTCT_transmission_rate)
        .def_readonly("PMTCT_dropout", &leapfrog::Children<double>::PMTCT_dropout)
        .def_readonly("PMTCT_input_is_percent", &leapfrog::Children<double>::PMTCT_input_is_percent)
        .def_readonly("breastfeeding_duration_art", &leapfrog::Children<double>::breastfeeding_duration_art)
        .def_readonly("breastfeeding_duration_no_art", &leapfrog::Children<double>::breastfeeding_duration_no_art)
        .def_readonly("mat_hiv_births", &leapfrog::Children<double>::mat_hiv_births)
        .def_readonly("mat_prev_input", &leapfrog::Children<double>::mat_prev_input)
        .def_readonly("prop_lt200", &leapfrog::Children<double>::prop_lt200)
        .def_readonly("prop_gte350", &leapfrog::Children<double>::prop_gte350)
        .def_readonly("incrate", &leapfrog::Children<double>::incrate)
        .def_readonly("ctx_val_is_percent", &leapfrog::Children<double>::ctx_val_is_percent)
        .def_readonly("hc_art_is_age_spec", &leapfrog::Children<double>::hc_art_is_age_spec)
        .def_readonly("hc_age_coarse", &leapfrog::Children<double>::hc_age_coarse)
        .def_readonly("abortion", &leapfrog::Children<double>::abortion)
        .def_readonly("patients_reallocated", &leapfrog::Children<double>::patients_reallocated)
        .def_readonly("hc_art_ltfu", &leapfrog::Children<double>::hc_art_ltfu)
        .def_readonly("hc_age_coarse_cd4", &leapfrog::Children<double>::hc_age_coarse_cd4)
        .def_readonly("adult_female_infections", &leapfrog::Children<double>::adult_female_infections)
        .def_readonly("adult_female_hivnpop", &leapfrog::Children<double>::adult_female_hivnpop)
        .def_readonly("total_births", &leapfrog::Children<double>::total_births);

    py::class_<leapfrog::Options<double>>(m, "Options")
        .def(py::init<int, int, int, int>())
        .def_readonly("hts_per_year", &leapfrog::Options<double>::hts_per_year)
        .def_readonly("dt", &leapfrog::Options<double>::dt)
        .def_readonly("p_idx_fertility_first", &leapfrog::Options<double>::p_idx_fertility_first)
        .def_readonly("p_fertility_age_groups", &leapfrog::Options<double>::p_fertility_age_groups)
        .def_readonly("p_idx_hiv_first_adult", &leapfrog::Options<double>::p_idx_hiv_first_adult)
        .def_readonly("adult_incidence_first_age_group", &leapfrog::Options<double>::adult_incidence_first_age_group)
        .def_readonly("pAG_INCIDPOP", &leapfrog::Options<double>::pAG_INCIDPOP)
        .def_readonly("ts_art_start", &leapfrog::Options<double>::ts_art_start)
        .def_readonly("hAG_15plus", &leapfrog::Options<double>::hAG_15plus)
        .def_readonly("hIDX_15PLUS", &leapfrog::Options<double>::hIDX_15PLUS);

    py::class_<leapfrog::DemographicProjectionParameters<double>>(m, "DemographicProjectionParameters")
        .def(py::init<const leapfrog::Demography<double>&>())
        .def_readonly("demography", &leapfrog::DemographicProjectionParameters<double>::demography);

    py::class_<leapfrog::HivSimulationParameters<double>>(m, "HivSimulationParameters")
        .def(py::init<const leapfrog::Incidence<double>&,
                      const leapfrog::NaturalHistory<double>&,
                      const leapfrog::Art<double>&>())
        .def_readonly("incidence", &leapfrog::HivSimulationParameters<double>::incidence)
        .def_readonly("natural_history", &leapfrog::HivSimulationParameters<double>::natural_history)
        .def_readonly("art", &leapfrog::HivSimulationParameters<double>::art);

    py::class_<leapfrog::ChildModelParameters<double>>(m, "ChildModelParameters")
        .def(py::init<const leapfrog::Children<double>&>())
        .def_readonly("children", &leapfrog::ChildModelParameters<double>::children);

    py::class_<leapfrog::Parameters<leapfrog::ChildModel, double>>(m, "Parameters")
        .def(py::init<const leapfrog::Options<double>&,
                      const leapfrog::DemographicProjectionParameters<double>&,
                      const leapfrog::HivSimulationParameters<double>&,
                      const leapfrog::ChildModelParameters<double>&>())
        .def_readonly("options", &leapfrog::Parameters<leapfrog::ChildModel, double>::options)
        .def_readonly("dp", &leapfrog::Parameters<leapfrog::ChildModel, double>::dp)
        .def_readonly("hiv", &leapfrog::Parameters<leapfrog::ChildModel, double>::hiv)
        .def_readonly("children", &leapfrog::Parameters<leapfrog::ChildModel, double>::children);

    py::class_<leapfrog::DemographicProjectionState<leapfrog::ChildModel, double, false>>(m, "DemographicProjectionState")
        .def(py::init<const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<1>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&>(),
                      py::arg("p_total_pop"),
                      py::arg("births"),
                      py::arg("p_total_pop_natural_deaths"))
        .def_readonly("p_total_pop", &leapfrog::DemographicProjectionState<leapfrog::ChildModel, double, false>::p_total_pop)
        .def_readonly("births", &leapfrog::DemographicProjectionState<leapfrog::ChildModel, double, false>::births)
        .def_readonly("p_total_pop_natural_deaths", &leapfrog::DemographicProjectionState<leapfrog::ChildModel, double, false>::p_total_pop_natural_deaths);

    py::class_<leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>>(m, "HivSimulationState")
        .def(py::init<const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hDS<leapfrog::ChildModel>, leapfrog::hAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hTS<leapfrog::ChildModel>, leapfrog::hDS<leapfrog::ChildModel>, leapfrog::hAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hDS<leapfrog::ChildModel>, leapfrog::hAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hTS<leapfrog::ChildModel>, leapfrog::hDS<leapfrog::ChildModel>, leapfrog::hAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hDS<leapfrog::ChildModel>, leapfrog::hAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&>(),
                      py::arg("p_hiv_pop"),
                      py::arg("p_hiv_pop_natural_deaths"),
                      py::arg("h_hiv_adult"),
                      py::arg("h_art_adult"),
                      py::arg("h_hiv_deaths_no_art"),
                      py::arg("p_infections"),
                      py::arg("h_hiv_deaths_art"),
                      py::arg("h_art_initiation"),
                      py::arg("p_hiv_deaths"))
        .def_readonly("p_hiv_pop", &leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>::p_hiv_pop)
        .def_readonly("p_hiv_pop_natural_deaths", &leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>::p_hiv_pop_natural_deaths)
        .def_readonly("h_hiv_adult", &leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>::h_hiv_adult)
        .def_readonly("h_art_adult", &leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>::h_art_adult)
        .def_readonly("h_hiv_deaths_no_art", &leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>::h_hiv_deaths_no_art)
        .def_readonly("p_infections", &leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>::p_infections)
        .def_readonly("h_hiv_deaths_art", &leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>::h_hiv_deaths_art)
        .def_readonly("h_art_initiation", &leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>::h_art_initiation)
        .def_readonly("p_hiv_deaths", &leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>::p_hiv_deaths);

    py::class_<leapfrog::ChildModelState<leapfrog::ChildModel, double, false>>(m, "ChildModelState")
        .def(py::init<const leapfrog::TensorType<double, Sizes<leapfrog::hc1DS<leapfrog::ChildModel>, leapfrog::hcTT<leapfrog::ChildModel>, leapfrog::hc1AG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hc2DS<leapfrog::ChildModel>, leapfrog::hcTT<leapfrog::ChildModel>, leapfrog::hc2AG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hTS<leapfrog::ChildModel>, leapfrog::hc1DS<leapfrog::ChildModel>, leapfrog::hc1AG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hTS<leapfrog::ChildModel>, leapfrog::hc2DS<leapfrog::ChildModel>, leapfrog::hc2AG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hc1DS<leapfrog::ChildModel>, leapfrog::hcTT<leapfrog::ChildModel>, leapfrog::hc1AG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hc2DS<leapfrog::ChildModel>, leapfrog::hcTT<leapfrog::ChildModel>, leapfrog::hc2AG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hTS<leapfrog::ChildModel>, leapfrog::hc1DS<leapfrog::ChildModel>, leapfrog::hc1AG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hTS<leapfrog::ChildModel>, leapfrog::hc2DS<leapfrog::ChildModel>, leapfrog::hc2AG<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<1>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hcAG_coarse<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hc1DS<leapfrog::ChildModel>, leapfrog::hcTT<leapfrog::ChildModel>, leapfrog::hcAG_end<leapfrog::ChildModel>, leapfrog::NS<leapfrog::ChildModel>>, false>&,
                      const leapfrog::TensorType<double, Sizes<1>, false>&,
                      const leapfrog::TensorType<double, Sizes<1>, false>&>(),
                      py::arg("hc1_hiv_pop"),
                      py::arg("hc2_hiv_pop"),
                      py::arg("hc1_art_pop"),
                      py::arg("hc2_art_pop"),
                      py::arg("hc1_noart_aids_deaths"),
                      py::arg("hc2_noart_aids_deaths"),
                      py::arg("hc1_art_aids_deaths"),
                      py::arg("hc2_art_aids_deaths"),
                      py::arg("hiv_births"),
                      py::arg("hc_art_init"),
                      py::arg("hc_art_need_init"),
                      py::arg("ctx_need"),
                      py::arg("ctx_mean"))
        .def_readonly("hc1_hiv_pop", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hc1_hiv_pop)
        .def_readonly("hc2_hiv_pop", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hc2_hiv_pop)
        .def_readonly("hc1_art_pop", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hc1_art_pop)
        .def_readonly("hc2_art_pop", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hc2_art_pop)
        .def_readonly("hc1_noart_aids_deaths", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hc1_noart_aids_deaths)
        .def_readonly("hc2_noart_aids_deaths", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hc2_noart_aids_deaths)
        .def_readonly("hc1_art_aids_deaths", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hc1_art_aids_deaths)
        .def_readonly("hc2_art_aids_deaths", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hc2_art_aids_deaths)
        .def_readonly("hiv_births", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hiv_births)
        .def_readonly("hc_art_init", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hc_art_init)
        .def_readonly("hc_art_need_init", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::hc_art_need_init)
        .def_readonly("ctx_need", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::ctx_need)
        .def_readonly("ctx_mean", &leapfrog::ChildModelState<leapfrog::ChildModel, double, false>::ctx_mean);

    py::class_<leapfrog::State<leapfrog::ChildModel, double, false>>(m, "State")
        .def(py::init<const leapfrog::DemographicProjectionState<leapfrog::ChildModel, double, false>&,
                      const leapfrog::HivSimulationState<leapfrog::ChildModel, double, false>&,
                      const leapfrog::ChildModelState<leapfrog::ChildModel, double, false>&>())
        .def_readonly("dp", &leapfrog::State<leapfrog::ChildModel, double, false>::dp)
        .def_readonly("hiv", &leapfrog::State<leapfrog::ChildModel, double, false>::hiv)
        .def_readonly("children", &leapfrog::State<leapfrog::ChildModel, double, false>::children);

    m.def(
        "project_single_year_cpp",
        &leapfrog::project_single_year<leapfrog::ChildModel, double, false>,
        py::arg("time_step"),
        py::arg("pars"),
        py::arg("state_curr"),
        py::arg("state_next"),
        "Project a single year of the model"
    );

    m.def(
        "set_initial_state_cpp",
        &leapfrog::set_initial_state<leapfrog::ChildModel, double, false>,
        py::arg("state"),
        py::arg("pars"),
        "Set initial state from the parameters"
    );

     m.def(
         "run_model_cpp",
         &leapfrog::simulate_model<leapfrog::ChildModel, double>,
         py::arg("params"),
         py::arg("proj_years"),
         py::arg("save_steps"),
         "Run a simulation model over a specified number of time steps"
     );
}
