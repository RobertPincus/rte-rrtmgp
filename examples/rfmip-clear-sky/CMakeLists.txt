#
# Example-specific utilities. Need to specify include directories because
# mo_rfmip_io.F90 uses netcdf directly
#
add_library(rfmip_clear_utils STATIC mo_rfmip_io.F90)
target_include_directories(
  rfmip_clear_utils
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)
target_link_libraries(rfmip_clear_utils PRIVATE shared-utils)

# ##############################################################################
#
# RRTMGP tests
#
# Build the two executables
#
foreach(
  test_executable IN
  ITEMS # cmake-format: sort
        rrtmgp_rfmip_lw
        rrtmgp_rfmip_sw
)
  add_executable(${test_executable} ${test_executable}.F90)
  target_include_directories(
    ${test_executable}
    PUBLIC
      $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
  )
  target_link_libraries(
    ${test_executable}
    PRIVATE rfmip_clear_utils
            rrtmgp-load-data
            shared-utils
            rte-rrtmgp::rte
            rte-rrtmgp::rrtmgp
            rte-rrtmgp::ssm
  )
endforeach()

#
# Copy empty output files
#
set(inoutputs
    rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
)

list(
  TRANSFORM inoutputs
  PREPEND ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/
          OUTPUT_VARIABLE inputs
)
#
# The tests write to the input files, therefore we copy them:
#
add_test(
  NAME copy_rrtmgp_rfmip_inputs
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different ${inputs}
    ${CMAKE_CURRENT_BINARY_DIR}/
)
set_tests_properties(
  copy_rrtmgp_rfmip_inputs
  PROPERTIES FIXTURES_REQUIRED
             fetch_rrtmgp_data
             FIXTURES_SETUP
             copy_rrtmgp_rfmip_inputs
)

#
# Comparision to RRTMGP RFMIP reference results
#
add_test(
  NAME run_rrtmgp_rfmip_lw
  COMMAND
    rrtmgp_rfmip_lw 8
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
)
set_tests_properties(
  run_rrtmgp_rfmip_lw
  PROPERTIES FIXTURES_REQUIRED
             "fetch_rrtmgp_data;copy_rrtmgp_rfmip_inputs"
             FIXTURES_SETUP
             run_rrtmgp_rfmip
)

add_test(
  NAME run_rrtmgp_rfmip_sw
  COMMAND
    rrtmgp_rfmip_sw 8
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
)
set_tests_properties(
  run_rrtmgp_rfmip_sw
  PROPERTIES FIXTURES_REQUIRED
             "fetch_rrtmgp_data;copy_rrtmgp_rfmip_inputs"
             FIXTURES_SETUP
             run_rrtmgp_rfmip
)

add_test(
  NAME check_rfmip_lw_sw
  COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/examples/compare-to-reference.py
    --ref_dir ${RRTMGP_DATA}/examples/rfmip-clear-sky/reference --tst_dir
    ${CMAKE_CURRENT_BINARY_DIR} --variables rld rlu rsd rsu --file_names
    ${inoutputs} --failure_threshold ${FAILURE_THRESHOLD}
)
set_tests_properties(
  check_rfmip_lw_sw
  PROPERTIES FIXTURES_REQUIRED "fetch_rrtmgp_data;run_rrtmgp_rfmip"
)
