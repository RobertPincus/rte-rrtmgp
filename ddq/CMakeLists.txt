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
# ##############################################################################
#
# Data
#
if(RTE_BUILD_TESTING)
  if(DDQ_DATA)
    add_test(NAME fetch_ddq_data COMMAND ${CMAKE_COMMAND} -E true)
    message(
      NOTICE
      "Using an external dataset from ${DDQ_DATA}: the data files will not be installed"
    )
  else()
    set(DDQ_DATA "${PROJECT_BINARY_DIR}/ddq-data")

    include(ExternalProject)
    ExternalProject_Add(
      ddq-data
      GIT_REPOSITORY https://github.com/ravikiran-hegde/ddq-data.git
      GIT_TAG "main"
      GIT_SHALLOW True
      EXCLUDE_FROM_ALL True
      PREFIX ddq-data-cmake
      SOURCE_DIR ${DDQ_DATA}
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
    )

    set(fetch_ddq_data_command
        ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR} --config
        "$<CONFIG>" --target ddq-data
    )

    add_test(NAME fetch_ddq_data COMMAND ${fetch_ddq_data_command})

    install(CODE "execute_process(COMMAND ${fetch_ddq_data_command})")
    install(
      FILES # cmake-format: sort
            ${DDQ_DATA}/gas-optics-lw.nc
            ${DDQ_DATA}/gas-optics-sw.nc
            ${DDQ_DATA}/rrtmgp-clouds-lw-ddq.nc
            ${DDQ_DATA}/rrtmgp-clouds-sw-ddq.nc
      DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/rte-rrtmgp/
    )

    set(REF_DIR examples/rte-examples/reference)

  endif()
endif()
