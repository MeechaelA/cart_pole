cmake_minimum_required(VERSION 3.16)

project(
  CartPole
  VERSION 1.0
  LANGUAGES CXX)

find_package(Eigen3 3.4 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)


include_directories(
  include 
  ${MPI_INCLUDE_PATH}
)

add_library(simulation_lib source/simulation.cpp)

add_library(cart_pole_lib source/cart_pole.cpp)

add_library(riccati_lib source/riccati_solver.cpp)
target_link_libraries(riccati_lib Eigen3::Eigen ${MPI_CXX_LIBRARIES})

target_link_libraries(simulation_lib cart_pole_lib riccati_lib Eigen3::Eigen OpenMP::OpenMP_CXX ${MPI_CXX_LIBRARIES})


add_executable(cart_pole source/main.cpp)

target_link_libraries(cart_pole PRIVATE simulation_lib Eigen3::Eigen ${MPI_CXX_LIBRARIES})