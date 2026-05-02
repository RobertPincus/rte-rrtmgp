add_library(
  fjx # cmake-format: sort
  mo_rte_fjx_interpolation_tables.F90
  mo_gas_optics_fjx.F90
  mo_cloud_optics_fjx.F90
  mo_fluxes_fjx.F90
)

add_library(rte-rrtmgp::fjx ALIAS fjx)

set_target_properties(
  fjx
  PROPERTIES EXPORT_NAME rte-rrtmgp::fjx
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
             INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
)

target_include_directories(
  fjx
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)

target_link_libraries(fjx PUBLIC rte-rrtmgp::rte)

install(TARGETS fjx EXPORT rte-rrtmgp-targets)
