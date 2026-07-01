# MCMC diagnostics tests

add_executable(mcmc_diagnostics_test mcmc_diagnostics_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME test_psrf COMMAND mcmc_diagnostics_test -tc=psrf)
add_test(NAME test_univariate_psrf COMMAND mcmc_diagnostics_test -tc=univariate_psrf)
add_test(NAME test_interval_psrf COMMAND mcmc_diagnostics_test -tc=interval_psrf)
add_test(NAME test_ess COMMAND mcmc_diagnostics_test -tc=ess)
add_test(NAME test_geweke COMMAND mcmc_diagnostics_test -tc=geweke)
add_test(NAME test_raftery COMMAND mcmc_diagnostics_test -tc=raftery)
TARGET_LINK_LIBRARIES(mcmc_diagnostics_test lp_solve ${MKL_LINK} coverage_config)
