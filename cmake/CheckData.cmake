# Verify that a set of data files is present in a directory.
#
# Run as: cmake -D DATA_DIR=<dir> \ -D DATA_FILES=<file1>[;<file2>...] \ -D
# DATA_LABEL=<label> \ -P CheckData.cmake
#
# Data sets are fetched at build time by their ExternalProject targets (part of
# the default ALL target). This script is used as a test fixture so that tests
# requiring the data fail with an actionable message if the project was not
# built before running ctest.

set(missing_files "")
foreach(data_file IN LISTS DATA_FILES)
  if(NOT EXISTS "${DATA_DIR}/${data_file}")
    list(APPEND missing_files "${DATA_DIR}/${data_file}")
  endif()
endforeach()

if(missing_files)
  string(REPLACE ";" "\n  " missing_pretty "${missing_files}")
  message(
    FATAL_ERROR
      "${DATA_LABEL} data files are missing:\n  ${missing_pretty}\n"
      "Build the project first (e.g. 'cmake --build <build-dir>') so the data "
      "is fetched before running ctest."
  )
endif()
