add_subdirectory(kernels)
add_subdirectory(frontend)

if(RRTMGP_DATA)
  add_test(NAME fetch_rrtmgp_data COMMAND ${CMAKE_COMMAND} -E true)
  message(
    NOTICE
    "Using an external dataset from ${RRTMGP_DATA}: the data files will not be installed"
  )
else()
  set(RRTMGP_DATA "${PROJECT_BINARY_DIR}/rrtmgp-data")

  include(ExternalProject)
  # Fetched at build time (part of the default ALL target) so the data is
  # present before install and before ctest.
  ExternalProject_Add(
    rrtmgp-data
    URL https://github.com/earth-system-radiation/rrtmgp-data/archive/refs/tags/v1.9.1.tar.gz
    PREFIX rrtmgp-data-cmake
    SOURCE_DIR ${RRTMGP_DATA}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD TRUE
  )

  # The build already fetched the data (rrtmgp-data is part of the ALL target).
  # This fixture verifies the files are present before any test that needs them
  # runs, and fails with an actionable message if the project was not built
  # first.
  add_test(
    NAME fetch_rrtmgp_data
    COMMAND
      ${CMAKE_COMMAND} -D DATA_DIR=${RRTMGP_DATA} -D DATA_LABEL=RRTMGP -D
      "DATA_FILES=rrtmgp-aerosols-merra-lw.nc;rrtmgp-aerosols-merra-sw.nc;rrtmgp-clouds-lw-bnd.nc;rrtmgp-clouds-sw-bnd.nc;rrtmgp-gas-lw-g128.nc;rrtmgp-gas-lw-g256.nc;rrtmgp-gas-sw-g112.nc;rrtmgp-gas-sw-g224.nc"
      -P ${PROJECT_SOURCE_DIR}/cmake/CheckData.cmake
  )

  install(
    FILES # cmake-format: sort
          ${RRTMGP_DATA}/rrtmgp-aerosols-merra-lw.nc
          ${RRTMGP_DATA}/rrtmgp-aerosols-merra-sw.nc
          ${RRTMGP_DATA}/rrtmgp-clouds-lw-bnd.nc
          ${RRTMGP_DATA}/rrtmgp-clouds-sw-bnd.nc
          ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
          ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
          ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
          ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/rte-rrtmgp/
  )
endif()

set_tests_properties(
  fetch_rrtmgp_data PROPERTIES FIXTURES_SETUP fetch_rrtmgp_data
)

# propagate the RRTMGP_DATA variable to the parent scope so that the sibling
# subdirectories (examples, tests) can inherit it
set(RRTMGP_DATA
    "${RRTMGP_DATA}"
    PARENT_SCOPE
)
