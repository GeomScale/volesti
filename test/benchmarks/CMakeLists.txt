# Benchmark executables (not run via ctest by default)

add_executable(benchmarks_sob benchmarks_sob.cpp)
add_executable(benchmarks_cg benchmarks_cg.cpp)
add_executable(benchmarks_cb benchmarks_cb.cpp)

TARGET_LINK_LIBRARIES(benchmarks_sob lp_solve ${MKL_LINK} coverage_config)
TARGET_LINK_LIBRARIES(benchmarks_cg lp_solve ${MKL_LINK} coverage_config)
TARGET_LINK_LIBRARIES(benchmarks_cb lp_solve ${MKL_LINK} coverage_config)

# CRHMC benchmarks (disabled by default)
# add_executable(benchmarks_crhmc benchmarks_crhmc.cpp)
# add_executable(benchmarks_crhmc_sampling benchmarks_crhmc_sampling.cpp)
# TARGET_LINK_LIBRARIES(benchmarks_crhmc_sampling lp_solve ${MKL_LINK} QD_LIB coverage_config)
# TARGET_LINK_LIBRARIES(benchmarks_crhmc lp_solve ${MKL_LINK} QD_LIB coverage_config)
