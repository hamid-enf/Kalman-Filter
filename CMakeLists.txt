cmake_minimum_required(VERSION 3.10)
project(kalman C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------

add_library(kalman STATIC
    kalman/src/kalman_matrix.c
    kalman/src/kalman_kf.c
    kalman/src/kalman_ekf.c
    kalman/src/kalman_ukf.c
    kalman/src/kalman_common.c
)

target_include_directories(kalman PUBLIC kalman/include)
target_compile_options(kalman PRIVATE -Wall -Wextra)

# Link the C math library where needed (not required on all platforms).
if(UNIX)
    target_link_libraries(kalman PUBLIC m)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

enable_testing()

add_executable(kalman_tests
    tests/test_framework.c
    tests/test_matrix.c
    tests/test_kf.c
    tests/test_ekf.c
    tests/test_ukf.c
    tests/test_extensions.c
    tests/main.c
)
target_include_directories(kalman_tests PRIVATE tests)
target_link_libraries(kalman_tests PRIVATE kalman)
add_test(NAME kalman_tests COMMAND kalman_tests)

# ---------------------------------------------------------------------------
# Examples (host)
# ---------------------------------------------------------------------------

set(KALMAN_EXAMPLES
    01_basic_1d
    02_position
    03_position_velocity
    04_position_velocity_acceleration
    05_multiple_measurements
    06_two_sensors
    07_predict_update_rates
    08_variable_dt
    09_ekf_basic
    10_ukf_basic
    11_outlier_gating
    12_adaptive_r
    13_rts_smoother
)

foreach(ex IN LISTS KALMAN_EXAMPLES)
    add_executable(example_${ex} examples/${ex}/main.c)
    target_include_directories(example_${ex} PRIVATE examples/common)
    target_link_libraries(example_${ex} PRIVATE kalman)
endforeach()

# ---------------------------------------------------------------------------
# Benchmark
# ---------------------------------------------------------------------------

add_executable(kalman_benchmark benchmarks/benchmark.c)
target_link_libraries(kalman_benchmark PRIVATE kalman)
