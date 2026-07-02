set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

add_library(
  ddq-load-data STATIC # cmake-format: sort
  mo_optics_utils_ddq.F90
)
target_link_libraries(
  ddq-load-data PUBLIC shared-utils rte-rrtmgp::ddq rte-rrtmgp::rte
)
target_include_directories(
  ddq-load-data
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
)
