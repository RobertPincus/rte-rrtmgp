cmake_minimum_required(VERSION 3.18)

project(
  rte-rrtmgp
  VERSION 1.9.2
  LANGUAGES Fortran
)

option(RTE_BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(RTE_BUILD_C_HEADERS "Build c headers" OFF)
option(RTE_BUILD_TESTING "Build tests" OFF)
option(RTE_ENABLE_SP "Enable single-precision floating-point model" OFF)

if(RTE_BUILD_C_HEADERS)
  set(RTE_BOOL "_Bool")
  set(RTE_PRECISION "double")
  set(RTE_EPSILON "2.2204460492503131e-16")

  if(RTE_ENABLE_SP)
    set(RTE_PRECISION "float")
    set(RTE_EPSILON "1.19209290E-07F")
  endif()
endif()

set(PREFERRED_KERNEL_MODES "default" "accel" "extern")
set(RTE_KERNEL_MODE
    "default"
    CACHE STRING "Select the RTE kernel mode: ${PREFERRED_KERNEL_MODES}"
)
set_property(CACHE RTE_KERNEL_MODE PROPERTY STRINGS ${PREFERRED_KERNEL_MODES})

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# CMake does not support the NEC compiler and usually mistakes it for GNU:
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # Check whether the Fortran compiler in use is NEC:
  include(CheckCompilerIsNEC)
  check_compiler_is_nec(RTE_RRTMGP_Fortran_is_NEC Fortran)
endif()

# GNUInstallDirs issues a warning if CMAKE_SIZEOF_VOID_P is not defined, which
# is the case with NAG and NEC. One way to circumvent that is to enable C
# language for the project:
if(CMAKE_Fortran_COMPILER_ID STREQUAL NAG OR RTE_RRTMGP_Fortran_is_NEC)
  enable_language(C)
endif()
include(GNUInstallDirs)

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU AND NOT RTE_RRTMGP_Fortran_is_NEC)
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

if(RTE_BUILD_TESTING OR RTE_BUILD_C_HEADERS)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

add_subdirectory(rte)
add_subdirectory(rrtmgp)
add_subdirectory(ssm)
add_subdirectory(ddq)

include(CTest)
if(RTE_BUILD_TESTING)
  if(NOT _RTE_RRTMGP_HAVE_PY_PACKAGES)
    include(CheckPython3Package)
    check_python3_package(numpy)
    check_python3_package(
      "netCDF4 or h5netcdf+scipy"
      CODE "try:
    import netCDF4
except:
    import h5netcdf
    import scipy
"
    )
    check_python3_package(
      "xarray>=0.12.2"
      CODE "import xarray
exit(tuple(map(int, xarray.__version__.split('.'))) < (0, 12, 2))"
    )
    check_python3_package(dask.array)
    set(_RTE_RRTMGP_HAVE_PY_PACKAGES
        TRUE
        CACHE INTERNAL
              "RTE-RRTMGP found all Python packages required for testing"
    )
  endif()

  find_package(NetCDF_Fortran REQUIRED)

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
        ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR} --config
        "$<CONFIG>" --target rrtmgp-data
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
  endif()

  set_tests_properties(
    fetch_rrtmgp_data PROPERTIES FIXTURES_SETUP fetch_rrtmgp_data
  )
  add_subdirectory(rrtmgp/data-loading-examples)
  add_subdirectory(examples)
  add_subdirectory(tests)
else()
  # Allow for 'make test' even if the tests are disabled:
  enable_testing()
endif()

export(
  EXPORT rte-rrtmgp-targets FILE ${PROJECT_BINARY_DIR}/rte-rrtmgp-targets.cmake
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/config.cmake.in
  ${PROJECT_BINARY_DIR}/rte-rrtmgp-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/rte-rrtmgp/cmake
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/rte-rrtmgp-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMinorVersion
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ TYPE INCLUDE)

install(
  EXPORT rte-rrtmgp-targets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/rte-rrtmgp/cmake
)

install(
  FILES ${PROJECT_BINARY_DIR}/rte-rrtmgp-config.cmake
        ${PROJECT_BINARY_DIR}/rte-rrtmgp-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/rte-rrtmgp/cmake
)
