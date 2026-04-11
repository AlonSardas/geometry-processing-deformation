add_executable(run_tests
  test_biharmonic_precompute.cpp
)

target_link_libraries(run_tests
  core
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(run_tests)
