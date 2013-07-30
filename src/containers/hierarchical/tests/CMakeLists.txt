set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Unit Tests
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Hierarchical Container:
add_executable(hierarchical_container_test container_test.cpp)
target_link_libraries(hierarchical_container_test gtest gtest_main)
add_test(hierarchical_container_test hierarchical_container_test)
