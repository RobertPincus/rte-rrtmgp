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
  # Fetched at build time (part of the default ALL target) so the data is
  # present before install and before ctest.
  ExternalProject_Add(
    rte-examples-data
    URL https://github.com/earth-system-radiation/rte-examples/archive/refs/tags/v0.1.1.tar.gz
    PREFIX rte-examples-data-cmake
    SOURCE_DIR ${RTE_EXAMPLES_DATA}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD TRUE
  )

  # The build already fetched the data (rte-examples-data is part of the ALL
  # target). This fixture verifies the files are present before any test that
  # needs them runs, and fails with an actionable message if the project was not
  # built first.
  add_test(
    NAME fetch_rte_examples_data
    COMMAND
      ${CMAKE_COMMAND} -D DATA_DIR=${RTE_EXAMPLES_DATA} -D
      DATA_LABEL=rte-examples -D
      "DATA_FILES=ckdmip-states.nc;rce-states.nc;rfmip-states.nc" -P
      ${PROJECT_SOURCE_DIR}/cmake/CheckData.cmake
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
          ddq-load-data
          shared-utils
          rte-rrtmgp::rte
          rte-rrtmgp::rrtmgp
          rte-rrtmgp::ssm
)

# ##############################################################################
#
# RTE tests: 2 schemes x 2 bands x 3 states = 12 tests
#
foreach(scheme IN ITEMS rrtmgp ssm ddq)
  if(scheme STREQUAL "rrtmgp")
    set(data_dir ${RRTMGP_DATA})
  elseif(scheme STREQUAL "ddq")
    set(data_dir ${DDQ_DATA})
  elseif(scheme STREQUAL "ssm")
    set(data_dir ${SSM_DATA})
  endif()

  foreach(band IN ITEMS lw sw)
    # scheme_data depends on scheme and band
    if(scheme STREQUAL "rrtmgp")
      if(band STREQUAL "lw")
        set(scheme_data ${data_dir}/rrtmgp-gas-lw-g128.nc)
      else()
        set(scheme_data ${data_dir}/rrtmgp-gas-sw-g112.nc)
      endif()
    elseif(scheme STREQUAL "ddq")
      if(band STREQUAL "lw")
        set(scheme_data ${data_dir}/gas_optics_lw.nc)
      else()
        set(scheme_data ${data_dir}/gas_optics_sw.nc)
      endif()
    elseif(scheme STREQUAL "ssm")
      set(scheme_data ${band})
    endif()

    foreach(state IN ITEMS rfmip ckdmip rce)
      set(output_data ${CMAKE_CURRENT_BINARY_DIR}/${scheme}-${band}-${state}.nc)
      add_test(
        NAME run_rte_examples_${scheme}_${band}_${state}
        COMMAND
          rte_examples 8 ${scheme} ${scheme_data}
          ${RTE_EXAMPLES_DATA}/${state}-states.nc ${output_data}
      )
      set_tests_properties(
        run_rte_examples_${scheme}_${band}_${state}
        PROPERTIES FIXTURES_REQUIRED
                   "fetch_${scheme}_data;fetch_rte_examples_data"
                   FIXTURES_SETUP
                   "run_rte_examples_${scheme}_${band}_${state}"
      )
    endforeach()
  endforeach()

  foreach(state IN ITEMS rfmip ckdmip rce)
    add_test(
      NAME check_rte_examples_${scheme}_${state}
      COMMAND
        ${Python3_EXECUTABLE}
        ${CMAKE_SOURCE_DIR}/examples/compare-to-reference.py --ref_dir
        ${data_dir}/examples/rte-examples/reference --tst_dir
        ${CMAKE_CURRENT_BINARY_DIR} --file_names
        ${CMAKE_CURRENT_BINARY_DIR}/${scheme}-lw-${state}.nc
        ${CMAKE_CURRENT_BINARY_DIR}/${scheme}-sw-${state}.nc
        --failure_threshold ${FAILURE_THRESHOLD}
    )
    set_tests_properties(
      check_rte_examples_${scheme}_${state}
      PROPERTIES
        FIXTURES_REQUIRED
        "fetch_ssm_data;run_rte_examples_${scheme}_lw_${state};run_rte_examples_${scheme}_sw_${state}"
    )
  endforeach()
endforeach()
