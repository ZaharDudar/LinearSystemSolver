project(TimeTests)
add_executable(TimeTestMulbyVector ./MulByVector.cpp)
target_link_libraries(TimeTestMulbyVector PRIVATE LinearSystemSolver)

add_executable(TimeTestIterMethods ./IterMethods.cpp)
target_link_libraries(TimeTestIterMethods PRIVATE LinearSystemSolver)

add_executable(TimeTestSpeededMPI ./SpeededMPI&others.cpp)
target_link_libraries(TimeTestSpeededMPI PRIVATE LinearSystemSolver)

add_executable(TestSORagainstOmega ./SORvsOmega.cpp)
target_link_libraries(TestSORagainstOmega PRIVATE LinearSystemSolver)

add_executable(TestUniversalBoost ./UniversalBoost.cpp)
target_link_libraries(TestUniversalBoost PRIVATE LinearSystemSolver)