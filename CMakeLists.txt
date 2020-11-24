cmake_minimum_required (VERSION 3.15)
project(nbody)

add_executable(${PROJECT_NAME} nbody.cpp)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	target_compile_options(${PROJECT_NAME} PRIVATE
		-fno-math-errno # sqrt prevents vectorization otherwise
		-mavx512f # if your system does not support his, use -march=native
	)
endif()
