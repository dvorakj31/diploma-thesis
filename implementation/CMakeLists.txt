cmake_minimum_required(VERSION 3.19)
project(dip)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic -fopenmp -lboost_mpi -lboost_serialization")


find_package(Boost COMPONENTS program_options serialization REQUIRED)

add_executable(dip src/main.cpp src/AbstractModel.h src/Options.h src/WeierstrassModel.cpp src/Lenstra.cpp src/EdwardsModel.cpp)
target_link_libraries(dip gmp ntl Boost::program_options Boost::serialization)


# Add tests
include(CTest)
add_subdirectory(tests)
