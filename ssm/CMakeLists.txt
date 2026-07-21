add_library(
  ssm # cmake-format: sort
  mo_optics_ssm.F90 mo_optics_ssm_kernels.F90
)

target_link_libraries(ssm PUBLIC rte)

target_include_directories(
  ssm
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
)

if(RTE_BUILD_SHARED_LIBS)
  set_target_properties(ssm PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

add_library(rte-rrtmgp::ssm ALIAS ssm)

set_target_properties(
  ssm
  PROPERTIES EXPORT_NAME rte-rrtmgp::ssm
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

install(TARGETS ssm EXPORT rte-rrtmgp-targets)

# ##############################################################################
#
# Fetch the rte-examples atmospheric state data
#
if(SSM_DATA)
  add_test(NAME fetch_ssm_data COMMAND ${CMAKE_COMMAND} -E true)
  message(
    NOTICE
    "Using an external dataset from ${SSM_DATA}: the data files will not be installed"
  )
else()
  set(SSM_DATA "${PROJECT_BINARY_DIR}/ssm-data")

  include(ExternalProject)
  # Fetched at build time (part of the default ALL target) so the data is
  # present before install and before ctest.
  ExternalProject_Add(
    ssm-data
    URL https://github.com/earth-system-radiation/ssm-data/archive/refs/tags/v0.0.1.tar.gz
    PREFIX ssm-data-cmake
    SOURCE_DIR ${SSM_DATA}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD TRUE
  )

  set(REF_DIR examples/rte-examples/reference)

  # The build already fetched the data (ssm-data is part of the ALL target).
  # This fixture verifies the files are present before any test that needs them
  # runs, and fails with an actionable message if the project was not built
  # first.
  add_test(
    NAME fetch_ssm_data
    COMMAND
      ${CMAKE_COMMAND} -D DATA_DIR=${SSM_DATA}/${REF_DIR} -D DATA_LABEL=SSM -D
      "DATA_FILES=ssm-lw-ckdmip.nc;ssm-lw-rce.nc;ssm-lw-rfmip.nc;ssm-sw-ckdmip.nc;ssm-sw-rce.nc;ssm-sw-rfmip.nc"
      -P ${PROJECT_SOURCE_DIR}/cmake/CheckData.cmake
  )
endif()

set_tests_properties(fetch_ssm_data PROPERTIES FIXTURES_SETUP fetch_ssm_data)

# propagate the SSM_DATA variable to the parent scope so that the sibling
# subdirectories (examples, tests) can inherit it
set(SSM_DATA
    "${SSM_DATA}"
    PARENT_SCOPE
)
