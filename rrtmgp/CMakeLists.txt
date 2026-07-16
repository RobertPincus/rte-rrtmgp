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
  ExternalProject_Add(
    rrtmgp-data
    GIT_REPOSITORY https://github.com/earth-system-radiation/rrtmgp-data.git
    GIT_TAG "v1.9"
    GIT_SHALLOW True
    EXCLUDE_FROM_ALL True
    PREFIX rrtmgp-data-cmake
    SOURCE_DIR ${RRTMGP_DATA}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )

  set(fetch_rrtmgp_data_command
      ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR} --config "$<CONFIG>"
      --target rrtmgp-data
  )

  add_test(NAME fetch_rrtmgp_data COMMAND ${fetch_rrtmgp_data_command})

  install(CODE "execute_process(COMMAND ${fetch_rrtmgp_data_command})")
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

  set_tests_properties(
    fetch_rrtmgp_data PROPERTIES FIXTURES_SETUP fetch_rrtmgp_data
  )
endif()
