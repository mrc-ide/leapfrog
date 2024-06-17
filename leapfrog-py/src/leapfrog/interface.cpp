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


    py::class_<leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>>(m, "BaseModelState")
        .def(py::init<const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&,
                      const leapfrog::TensorType<double, Sizes<1>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hDS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hTS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hDS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hDS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hTS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hDS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::hDS<leapfrog::BaseModelFullAgeStratification>, leapfrog::hAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&,
                      const leapfrog::TensorType<double, Sizes<leapfrog::pAG<leapfrog::BaseModelFullAgeStratification>, leapfrog::NS<leapfrog::BaseModelFullAgeStratification>>, false>&>(),
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
        .def("reset", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::reset)
        .def_readonly("p_total_pop", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::p_total_pop)
        .def_readonly("births", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::births)
        .def_readonly("p_total_pop_natural_deaths", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::p_total_pop_natural_deaths)
        .def_readonly("p_hiv_pop", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::p_hiv_pop)
        .def_readonly("p_hiv_pop_natural_deaths", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::p_hiv_pop_natural_deaths)
        .def_readonly("h_hiv_adult", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::h_hiv_adult)
        .def_readonly("h_art_adult", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::h_art_adult)
        .def_readonly("h_hiv_deaths_no_art", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::h_hiv_deaths_no_art)
        .def_readonly("p_infections", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::p_infections)
        .def_readonly("h_hiv_deaths_art", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::h_hiv_deaths_art)
        .def_readonly("h_art_initiation", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::h_art_initiation)
        .def_readonly("p_hiv_deaths", &leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>::p_hiv_deaths);

    py::class_<leapfrog::ChildModelState<leapfrog::BaseModelFullAgeStratification, double, false>>(m, "ChildModelState")
        .def(py::init<>())
        .def("reset", &leapfrog::ChildModelState<leapfrog::BaseModelFullAgeStratification, double, false>::reset);

    py::class_<leapfrog::State<leapfrog::BaseModelFullAgeStratification, double, false>>(m, "State")
        .def(py::init<const leapfrog::BaseModelState<leapfrog::BaseModelFullAgeStratification, double, false>&,
                      const leapfrog::ChildModelState<leapfrog::BaseModelFullAgeStratification, double, false>&>())
        .def("reset", &leapfrog::State<leapfrog::BaseModelFullAgeStratification, double, false>::reset)
        .def_readonly("base", &leapfrog::State<leapfrog::BaseModelFullAgeStratification, double, false>::base)
        .def_readonly("children", &leapfrog::State<leapfrog::BaseModelFullAgeStratification, double, false>::children);

    m.def(
        "project_single_year_cpp",
        &leapfrog::project_single_year<leapfrog::BaseModelFullAgeStratification, double, false>,
        py::arg("time_step"),
        py::arg("hiv_steps"),
        py::arg("parameters"),
        py::arg("state_curr"),
        py::arg("state_next"),
        "Project a single year of the model"
    );

    m.def(
        "set_initial_state_cpp",
        &leapfrog::set_initial<leapfrog::BaseModelFullAgeStratification, double, false>,
        py::arg("state"),
        py::arg("parameters"),
        "Set initial state from the parameters"
    );

     m.def(
         "run_model_cpp",
         &leapfrog::simulate_model<leapfrog::BaseModelFullAgeStratification, double>,
         py::arg("parameters"),
         py::arg("proj_years"),
         py::arg("hiv_steps"),
         py::arg("save_steps"),
         "Run a simulation model over a specified number of time steps"
     );
}
