cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

# Add Eigen
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

pybind11_add_module(leapfrog src/interface/interface.cpp)
target_include_directories(leapfrog PRIVATE ../inst/include)
target_link_libraries (leapfrog PRIVATE Eigen3::Eigen)

install(TARGETS leapfrog LIBRARY DESTINATION .)