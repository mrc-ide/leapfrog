cmake_minimum_required(VERSION 3.19)
project(standalone_model)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
add_executable(simulate_model simulate_model.cpp)
target_link_libraries(simulate_model HDF5::HDF5)
target_include_directories(simulate_model PRIVATE ../leapfrogr/inst/include/)

