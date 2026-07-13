#
# Example-specific utilities. Need to specify include directories because
# mo_rte_examples_io.F90 uses netcdf directly
#
set(executable rte_examples)

add_library(${executable}_io STATIC mo_${executable}_io.F90)
target_include_directories(
  ${executable}_io
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)
target_link_libraries(rte_examples_io PRIVATE shared-utils)

# ##############################################################################
#
# Fetch the rte-examples atmospheric state data
#
if(RTE_EXAMPLES_DATA)
  add_test(NAME fetch_rte_examples_data COMMAND ${CMAKE_COMMAND} -E true)
  message(
    NOTICE
    "Using an external dataset from ${RTE_EXAMPLES_DATA}: the data files will not be installed"
  )
else()
  set(RTE_EXAMPLES_DATA "${PROJECT_BINARY_DIR}/rte-examples-data")

  include(ExternalProject)
  ExternalProject_Add(
    rte-examples-data
    GIT_REPOSITORY https://github.com/earth-system-radiation/rte-examples.git
    GIT_TAG "v0.1.1"
    GIT_SHALLOW True
    EXCLUDE_FROM_ALL True
    PREFIX rte-examples-data-cmake
    SOURCE_DIR ${RTE_EXAMPLES_DATA}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )

  set(fetch_rte_examples_data_command
      ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --config "$<CONFIG>"
      --target rte-examples-data
  )

  add_test(
    NAME fetch_rte_examples_data COMMAND ${fetch_rte_examples_data_command}
  )

  install(CODE "execute_process(COMMAND ${fetch_rte_examples_data_command})")
  install(
    FILES # cmake-format: sort
          ${RTE_EXAMPLES_DATA}/ckdmip-states.nc
          ${RTE_EXAMPLES_DATA}/rce-states.nc
          ${RTE_EXAMPLES_DATA}/rfmip-states.nc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/rte-rrtmgp/
  )
endif()

set_tests_properties(
  fetch_rte_examples_data PROPERTIES FIXTURES_SETUP fetch_rte_examples_data
)

# ##############################################################################
#
# RTE tests using homogenized examples
#
# Build the two executables
#
add_executable(${executable} ${executable}.F90)
target_include_directories(
  ${executable}
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
)
target_link_libraries(
  ${executable}
  PRIVATE ${executable}_io
          rrtmgp-load-data
          shared-utils
          rte-rrtmgp::rte
          rte-rrtmgp::rrtmgp
          rte-rrtmgp::ssm
)

# ##############################################################################
#
# RTE tests: 2 schemes x 2 bands x 3 states = 12 tests
#
foreach(scheme IN ITEMS rrtmgp ssm)
  foreach(band IN ITEMS lw sw)
    # scheme_data depends on scheme and band
    if(scheme STREQUAL "rrtmgp")
      if(band STREQUAL "lw")
        set(scheme_data ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc)
      else()
        set(scheme_data ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc)
      endif()
    else() # ssm
      set(scheme_data ${band})
    endif()

    foreach(state IN ITEMS rfmip ckdmip rce)
      add_test(
        NAME run_rte_examples_${scheme}_${band}_${state}
        COMMAND
          rte_examples 8 ${scheme} ${scheme_data}
          ${RTE_EXAMPLES_DATA}/${state}-states.nc
          ${CMAKE_CURRENT_BINARY_DIR}/${scheme}_${band}_${state}.nc
      )
      set_tests_properties(
        run_rte_examples_${scheme}_${band}_${state}
        PROPERTIES FIXTURES_REQUIRED
                   "fetch_rrtmgp_data;fetch_rte_examples_data"
      )
    endforeach()
  endforeach()
endforeach()
