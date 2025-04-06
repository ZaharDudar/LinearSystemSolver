include(AddGoogleTest)

add_executable(
  unit_tests
  TestMatrix.cpp
  TestCSRMatrix.cpp
  TestSolving.cpp
  TestVectorOps.cpp
  TestIterMethods.cpp
)

target_link_libraries(
  unit_tests PRIVATE LinearSystemSolver GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(unit_tests)