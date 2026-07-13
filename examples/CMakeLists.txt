if(RTE_ENABLE_SP)
  set(default_FAILURE_THRESHOLD "3.5e-1")
else()
  set(default_FAILURE_THRESHOLD "7.e-4")
endif()

set(FAILURE_THRESHOLD
    "${default_FAILURE_THRESHOLD}"
    CACHE STRING "Default failure threshold for tests"
)

add_subdirectory(shared-utils)
add_subdirectory(all-sky)
add_subdirectory(rfmip-clear-sky)
add_subdirectory(rte-examples)
