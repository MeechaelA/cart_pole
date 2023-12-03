cmake_minimum_required(VERSION 3.16)

project(
  CartPole
  VERSION 1.0
  LANGUAGES CXX
)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP REQUIRED)


include_directories(
  include 
  ${MPI_INCLUDE_PATH}
)

include_directories(dependencies)
add_subdirectory(dependencies)

add_library(simulation_lib source/simulation.cpp)

add_library(cart_pole_lib source/cart_pole.cpp)

add_library(riccati_lib source/riccati_solver.cpp)
target_link_libraries(riccati_lib)

target_link_libraries(simulation_lib cart_pole_lib riccati_lib OpenMP::OpenMP_CXX)

add_executable(cart_pole source/main.cpp)

target_link_libraries(cart_pole PRIVATE simulation_lib nlohmann_json::nlohmann_json)