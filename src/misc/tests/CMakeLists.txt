set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Unit Tests
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Types
add_executable(integer integer_test.cpp)
target_link_libraries(integer gtest gtest_main)
add_test(integer integer)
