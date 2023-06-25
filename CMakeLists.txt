cmake_minimum_required(VERSION 3.10)

project(excitation_ratios Fortran)
cmake_policy(SET CMP0097 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_Fortran_FLAGS -fallow-argument-mismatch)

add_executable(${PROJECT_NAME} 
                src/main_ER.f
                src/forward_model_mineos.f
                src/minos_bran.f)
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

install(TARGETS ${PROJECT_NAME})