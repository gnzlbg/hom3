set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Unit Tests
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Unit Test: container.hpp
add_executable(container_test container_test.cpp)
target_link_libraries(container_test gtest gtest_main ${HOM3_BASE_LIBS})
add_test(container_test container_test)
