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
if(RTE_BUILD_TESTING)
  if(SSM_DATA)
    add_test(NAME fetch_ssm_data COMMAND ${CMAKE_COMMAND} -E true)
    message(
      NOTICE
      "Using an external dataset from ${SSM_DATA}: the data files will not be installed"
    )
  else()
    set(SSM_DATA "${PROJECT_BINARY_DIR}/ssm-data")

    include(ExternalProject)
    ExternalProject_Add(
      ssm-data
      GIT_REPOSITORY https://github.com/earth-system-radiation/ssm-data.git
      GIT_TAG "main"
      GIT_SHALLOW True
      EXCLUDE_FROM_ALL True
      PREFIX ssm-data-cmake
      SOURCE_DIR ${SSM_DATA}
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
    )

    set(fetch_ssm_data_command
        ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --config "$<CONFIG>"
        --target ssm-data
    )

    add_test(NAME fetch_ssm_data COMMAND ${fetch_ssm_data_command})

    set(REF_DIR examples/rte-examples/reference)

    install(CODE "execute_process(COMMAND ${fetch_ssm_data_command})")
    install(
      FILES # cmake-format: sort
            ${SSM_DATA}/${REF_DIR}/ssm-lw-ckdmip.nc
            ${SSM_DATA}/${REF_DIR}/ssm-lw-rce.nc
            ${SSM_DATA}/${REF_DIR}/ssm-lw-rfmip.nc
            ${SSM_DATA}/${REF_DIR}/ssm-sw-ckdmip.nc
            ${SSM_DATA}/${REF_DIR}/ssm-sw-rce.nc
            ${SSM_DATA}/${REF_DIR}/ssm-sw-rfmip.nc
      DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/rte-rrtmgp/
    )
  endif()

  set_tests_properties(fetch_ssm_data PROPERTIES FIXTURES_SETUP fetch_ssm_data)
endif()
