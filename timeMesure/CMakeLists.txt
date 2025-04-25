project(TimeTests)
add_executable(TimeTestMulbyVector ./MulByVector.cpp)
target_link_libraries(TimeTestMulbyVector PRIVATE LinearSystemSolver)

add_executable(TimeTestIterMethods ./IterMethods.cpp)
target_link_libraries(TimeTestIterMethods PRIVATE LinearSystemSolver)

add_executable(TimeTestSpeededMPI ./SpeededMPI&others.cpp)
target_link_libraries(TimeTestSpeededMPI PRIVATE LinearSystemSolver)
target_compile_options(TimeTestSpeededMPI PRIVATE -O3)

add_executable(TestSORagainstOmega ./SORvsOmega.cpp)
target_link_libraries(TestSORagainstOmega PRIVATE LinearSystemSolver)

add_executable(TestUniversalBoost ./UniversalBoost.cpp)
target_link_libraries(TestUniversalBoost PRIVATE LinearSystemSolver)
target_compile_options(TestUniversalBoost PRIVATE -O3)


add_executable(TestCG-BiCG-CSG_Sym ./CG-BiCG-CGS_Sym.cpp)
target_link_libraries(TestCG-BiCG-CSG_Sym PRIVATE LinearSystemSolver)
target_compile_options(TestCG-BiCG-CSG_Sym PRIVATE -O3)

add_executable(TestBiCG-CSG_NotSym ./BiCG-CGS_NotSym.cpp)
target_link_libraries(TestBiCG-CSG_NotSym PRIVATE LinearSystemSolver)
target_compile_options(TestBiCG-CSG_NotSym PRIVATE -O3)