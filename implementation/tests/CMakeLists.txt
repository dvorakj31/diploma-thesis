add_definitions(-DBOOST_ASIO_DISABLE_CONCEPTS)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

add_executable(
        test_lenstra
        test_lenstra.cpp
        test_runner.cpp
        ../src/AbstractModel.h
        ../src/Options.h
        ../src/WeierstrassModel.cpp
        ../src/Lenstra.cpp
        ../src/EdwardsModel.cpp
)

target_link_libraries(test_lenstra ntl Boost::unit_test_framework)

enable_testing()
add_test(test_lenstra test_lenstra)