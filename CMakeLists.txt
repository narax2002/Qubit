cmake_minimum_required(VERSION 3.16)
project(Qubit)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

add_compile_options(-fopenmp -Wall)

include_directories(${PROJECT_SOURCE_DIR}/include/)
link_directories(${PROJECT_SOURCE_DIR}/lib)

set(src
	src/qubit_operator.cpp
	src/one_qubit_gates.cpp
	src/multi_qubit_gates.cpp
	src/qubit_Bernstein_Vazirani.cpp
	src/qubit_QFT.cpp
	src/qubit_Shor.cpp
	)

add_library(Qubit SHARED ${src})

add_executable(Sample sample.cpp)
target_link_libraries(Sample Qubit gomp)