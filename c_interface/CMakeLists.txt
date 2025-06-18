cmake_minimum_required(VERSION 3.10)

project(LeapfrogDll)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the source files to the project
set(SOURCE_FILES
    c_interface.cpp
)

# Create the shared library (DLL)
add_library(leapfrog SHARED ${SOURCE_FILES})

message(CMAKE_CURRENT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

target_include_directories(leapfrog
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../r-package/inst/include
)

find_package(Eigen3 REQUIRED)
target_link_libraries(leapfrog Eigen3::Eigen)
