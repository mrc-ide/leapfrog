#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen/tensor.h>
#include "../../inst/include/types.hpp"
#include "../../inst/include/state_types.hpp"
#include "../../inst/include/leapfrog-py.hpp"
#include "../../inst/include/frogger.hpp"

namespace py = pybind11;

namespace {
using Eigen::Sizes;
using Eigen::TensorFixedSize;
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
                      const leapfrog::TensorMap1<double>&,
                      const leapfrog::TensorMap2<double>&,
                      const leapfrog::TensorMap2<int>&,
                      const leapfrog::Tensor1<double>&,
                      double>())
        .def_readonly("idx_hm_elig", &leapfrog::Art<double>::idx_hm_elig)
        .def_readonly("mortality", &leapfrog::Art<double>::mortality)
        .def_readonly("mortaility_time_rate_ratio", &leapfrog::Art<double>::mortaility_time_rate_ratio)
        .def_readonly("dropout", &leapfrog::Art<double>::dropout)
        .def_readonly("adults_on_art", &leapfrog::Art<double>::adults_on_art)
        .def_readonly("adults_on_art_is_percent", &leapfrog::Art<double>::adults_on_art_is_percent)
        .def_readonly("h_art_stage_dur", &leapfrog::Art<double>::h_art_stage_dur)
        .def_readonly("initiation_mortality_weight", &leapfrog::Art<double>::initiation_mortality_weight);

//    TODO: Add child model support
//    py::class_<Children<double>>(leapfrog, "Children")
//            .def(py::init<const TensorMap1<double>&,
//                      const TensorMap1<double>&,
//                      const TensorMap2<double>&,
//                      const TensorMap3<double>&,
//                      const TensorMap3<double>&,
//                      const TensorMap1<double>&,
//                      const TensorMap1<double>&,
//                      double,
//                      const TensorMap1<double>&,
//                      const TensorMap1<int>&,
//                      const Tensor2<double>&,
//                      const TensorMap3<double>&,
//                      const TensorMap3<double>&,
//                      const TensorMap3<double>&,
//                      const TensorMap1<int>&,
//                      const TensorMap2<double>&,
//                      const TensorMap2<double>&,
//                      const TensorMap2<double>&,
//                      const TensorMap1<double>&,
//                      const TensorMap1<double>&,
//                      const TensorMap1<double>&,
//                      const TensorMap1<double>&,
//                      double,
//                      const TensorMap2<double>&,
//                      const TensorMap2<double>&,
//                      const TensorMap3<double>&,
//                      const TensorMap2<double>&,
//                      const TensorMap1<int>&,
//                      const TensorMap2<double>&,
//                      const TensorMap2<double>&,
//                      const TensorMap1<double>&,
//                      const TensorMap1<int>&,
//                      const TensorMap1<double>&,
//                      const TensorMap1<double>&,
//                      const TensorMap1<int>&,
//                      const TensorMap1<int>&,
//                      const TensorMap1<double>&>())
//        .def_readonly("hc_nosocomial", &Children<double>::hc_nosocomial)
//        .def_readonly("hc1_cd4_dist", &Children<double>::hc1_cd4_dist)
//        .def_readonly("hc_cd4_transition", &Children<double>::hc_cd4_transition)
//        .def_readonly("hc1_cd4_mort", &Children<double>::hc1_cd4_mort)
//        .def_readonly("hc2_cd4_mort", &Children<double>::hc2_cd4_mort)
//        .def_readonly("hc1_cd4_prog", &Children<double>::hc1_cd4_prog)
//        .def_readonly("hc2_cd4_prog", &Children<double>::hc2_cd4_prog)
//        .def_readonly("ctx_effect", &Children<double>::ctx_effect)
//        .def_readonly("ctx_val", &Children<double>::ctx_val)
//        .def_readonly("hc_art_elig_age", &Children<double>::hc_art_elig_age)
//        .def_readonly("hc_art_elig_cd4", &Children<double>::hc_art_elig_cd4)
//        .def_readonly("hc_art_mort_rr", &Children<double>::hc_art_mort_rr)
//        .def_readonly("hc1_art_mort", &Children<double>::hc1_art_mort)
//        .def_readonly("hc2_art_mort", &Children<double>::hc2_art_mort)
//        .def_readonly("hc_art_isperc", &Children<double>::hc_art_isperc)
//        .def_readonly("hc_art_val", &Children<double>::hc_art_val)
//        .def_readonly("hc_art_init_dist", &Children<double>::hc_art_init_dist)
//        .def_readonly("adult_cd4_dist", &Children<double>::adult_cd4_dist)
//        .def_readonly("fert_mult_by_age", &Children<double>::fert_mult_by_age)
//        .def_readonly("fert_mult_off_art", &Children<double>::fert_mult_off_art)
//        .def_readonly("fert_mult_on_art", &Children<double>::fert_mult_on_art)
//        .def_readonly("total_fertility_rate", &Children<double>::total_fertility_rate)
//        .def_readonly("local_adj_factor", &Children<double>::local_adj_factor)
//        .def_readonly("PMTCT", &Children<double>::PMTCT)
//        .def_readonly("vertical_transmission_rate", &Children<double>::vertical_transmission_rate)
//        .def_readonly("PMTCT_transmission_rate", &Children<double>::PMTCT_transmission_rate)
//        .def_readonly("PMTCT_dropout", &Children<double>::PMTCT_dropout)
//        .def_readonly("PMTCT_input_is_percent", &Children<double>::PMTCT_input_is_percent)
//        .def_readonly("breastfeeding_duration_art", &Children<double>::breastfeeding_duration_art)
//        .def_readonly("breastfeeding_duration_no_art", &Children<double>::breastfeeding_duration_no_art)
//        .def_readonly("mat_hiv_births", &Children<double>::mat_hiv_births)
//        .def_readonly("mat_prev_input", &Children<double>::mat_prev_input)
//        .def_readonly("prop_lt200", &Children<double>::prop_lt200)
//        .def_readonly("prop_gte350", &Children<double>::prop_gte350)
//        .def_readonly("incrate", &Children<double>::incrate)
//        .def_readonly("ctx_val_is_percent", &Children<double>::ctx_val_is_percent)
//        .def_readonly("hc_art_is_age_spec", &Children<double>::hc_art_is_age_spec)
//        .def_readonly("hc_age_coarse", &Children<double>::hc_age_coarse);

    py::class_<leapfrog::Options<double>>(m, "Options")
        .def(py::init<int, int, int>())
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

    py::class_<leapfrog::BaseModelParameters<double>>(m, "BaseModelParameters")
        .def(py::init<const leapfrog::Options<double>&,
                      const leapfrog::Demography<double>&,
                      const leapfrog::Incidence<double>&,
                      const leapfrog::NaturalHistory<double>&,
                      const leapfrog::Art<double>&>())
        .def_readonly("options", &leapfrog::BaseModelParameters<double>::options)
        .def_readonly("demography", &leapfrog::BaseModelParameters<double>::demography)
        .def_readonly("incidence", &leapfrog::BaseModelParameters<double>::incidence)
        .def_readonly("natural_history", &leapfrog::BaseModelParameters<double>::natural_history)
        .def_readonly("art", &leapfrog::BaseModelParameters<double>::art);


//    py::class_<ChildModelParameters<ChildModel, double>>(m, "ChildModelParameters")
//        .def_readonly("children", &ChildModelParameters<ChildModel, double>::children);

    py::class_<leapfrog::ChildModelParameters<leapfrog::BaseModelFullAgeStratification, double>>(m, "BaseModelChildParameters")
        .def(py::init<>());

    py::class_<leapfrog::Parameters<leapfrog::BaseModelFullAgeStratification, double>>(m, "Parameters")
        .def(py::init<const leapfrog::BaseModelParameters<double>&,
                      const leapfrog::ChildModelParameters<leapfrog::BaseModelFullAgeStratification, double>&>())
        .def_readonly("base", &leapfrog::Parameters<leapfrog::BaseModelFullAgeStratification, double>::base)
        .def_readonly("children", &leapfrog::Parameters<leapfrog::BaseModelFullAgeStratification, double>::children);

    py::class_<leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>>(m, "BaseModelState")
        .def(py::init<const leapfrog::Parameters<leapfrog::BaseModelFullAgeStratification, double>&>())
        .def(py::init<const TensorFixedSize<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&,
                      double,
                      TensorFixedSize<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&,
                      TensorFixedSize<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&,
                      TensorFixedSize<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&,
                      TensorFixedSize<double, Sizes<leapfrog::hDS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&,
                      TensorFixedSize<double, Sizes<leapfrog::hTS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hDS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&,
                      TensorFixedSize<double, Sizes<leapfrog::hDS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&,
                      TensorFixedSize<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&,
                      TensorFixedSize<double, Sizes<leapfrog::hTS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hDS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&,
                      TensorFixedSize<double, Sizes<leapfrog::hDS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&,
                      TensorFixedSize<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>>&>(),
                      py::arg("p_total_pop"),
                      py::arg("births"),
                      py::arg("p_total_pop_natural_deaths"),
                      py::arg("p_hiv_pop"),
                      py::arg("p_hiv_pop_natural_deaths"),
                      py::arg("h_hiv_adult"),
                      py::arg("h_art_adult"),
                      py::arg("h_hiv_deaths_no_art"),
                      py::arg("p_infections"),
                      py::arg("h_hiv_deaths_art"),
                      py::arg("h_art_initiation"),
                      py::arg("p_hiv_deaths"))
        .def("reset", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::reset)
        .def_readonly("p_total_pop", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::p_total_pop)
        .def_readonly("births", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::births)
        .def_readonly("p_total_pop_natural_deaths", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::p_total_pop_natural_deaths)
        .def_readonly("p_hiv_pop", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::p_hiv_pop)
        .def_readonly("p_hiv_pop_natural_deaths", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::p_hiv_pop_natural_deaths)
        .def_readonly("h_hiv_adult", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::h_hiv_adult)
        .def_readonly("h_art_adult", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::h_art_adult)
        .def_readonly("h_hiv_deaths_no_art", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::h_hiv_deaths_no_art)
        .def_readonly("p_infections", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::p_infections)
        .def_readonly("h_hiv_deaths_art", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::h_hiv_deaths_art)
        .def_readonly("h_art_initiation", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::h_art_initiation)
        .def_readonly("p_hiv_deaths", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>::p_hiv_deaths);

    py::class_<leapfrog::ChildModelState<leapfrog::BaseModelFullAgeStratification, double>>(m, "ChildModelState")
        .def(py::init<const leapfrog::Parameters<leapfrog::BaseModelFullAgeStratification, double>&>())
        .def("reset", &leapfrog::ChildModelState<leapfrog::BaseModelFullAgeStratification, double>::reset);

    py::class_<leapfrog::State<leapfrog::BaseModelFullAgeStratification, double>>(m, "State")
        .def(py::init<const leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double>&,
                      const leapfrog::ChildModelState<leapfrog::BaseModelFullAgeStratification, double>&>())
        .def(py::init<const leapfrog::Parameters<leapfrog::BaseModelFullAgeStratification, double>&>())
        .def("reset", &leapfrog::State<leapfrog::BaseModelFullAgeStratification, double>::reset)
        .def_readonly("base", &leapfrog::State<leapfrog::BaseModelFullAgeStratification, double>::base)
        .def_readonly("children", &leapfrog::State<leapfrog::BaseModelFullAgeStratification, double>::children);

    m.def(
        "project_single_year_cpp",
        &leapfrog::project_single_year<leapfrog::BaseModelFullAgeStratification, double>,
        py::arg("time_step"),
        py::arg("pars"),
        py::arg("state_curr"),
        py::arg("state_next"),
        "Project a single year of the model"
    );

    m.def(
        "set_initial_state_cpp",
        &leapfrog::set_initial_state<leapfrog::BaseModelFullAgeStratification, double>,
        py::arg("state"),
        py::arg("pars"),
        "Set initial state from the parameters"
    );

    m.def(
        "run_model_cpp",
        &leapfrog::simulate_model<leapfrog::BaseModelFullAgeStratification, double>,
        py::arg("params"),
        py::arg("proj_years"),
        py::arg("save_steps"),
        "Run a simulation model over a specified number of time steps"
    );
}