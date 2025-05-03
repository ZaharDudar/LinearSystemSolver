include(AddGoogleTest)

add_executable(
  unit_tests
  TestMatrix.cpp
  TestColMajMatrix.cpp
  TestCSRMatrix.cpp
  TestSolving.cpp
  TestVectorOps.cpp
  TestIterMethods.cpp
  TestGradientMethods.cpp
  TestGMRES.cpp
)

target_link_libraries(
  unit_tests PRIVATE LinearSystemSolver GTest::gtest_main
)
target_compile_options(unit_tests PRIVATE -g)

include(GoogleTest)
gtest_discover_tests(unit_tests)