cmake_minimum_required(VERSION 3.16)

project(
  CartPole
  VERSION 1.0
  LANGUAGES CXX)

include_directories(include)

add_library(simulation_lib source/simulation.cpp)

add_library(cart_pole_lib source/cart_pole.cpp)
target_link_libraries(simulation_lib cart_pole_lib)


find_package(drake CONFIG REQUIRED PATHS /opt/drake)

get_filename_component(DRAKE_PYTHONPATH "${drake_DIR}" DIRECTORY)
get_filename_component(DRAKE_PYTHONPATH "${DRAKE_PYTHONPATH}" DIRECTORY)
set(DRAKE_PYTHONPATH
  "${DRAKE_PYTHONPATH}/python${PYTHON_VERSION_MAJOR_MINOR}/site-packages"
)


add_executable(cart_pole source/main.cpp)

target_link_libraries(cart_pole PRIVATE simulation_lib)