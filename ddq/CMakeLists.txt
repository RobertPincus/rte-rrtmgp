add_library(
  ddq # cmake-format: sort
  mo_gas_optics_ddq.F90 mo_gas_optics_ddq_kernels.F90
)

add_library(rte-rrtmgp::ddq ALIAS ddq)

set_target_properties(
  ddq
  PROPERTIES EXPORT_NAME rte-rrtmgp::ddq
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
             INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
)

target_include_directories(
  ddq
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
)

target_link_libraries(ddq PRIVATE rte)

install(TARGETS ddq EXPORT rte-rrtmgp-targets)
