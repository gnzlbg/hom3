set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Unit Tests
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Grid Node Container:
# add_executable(grid_container_test container_test.cpp)
# target_link_libraries(grid_container_test gtest gtest_main)
# add_test(grid_container_test grid_container_test)

# # Grid:
# add_executable(grid_test grid_test.cpp)
# target_link_libraries(grid_test gtest gtest_main)
# add_test(grid_test grid_test)

# Debug grid
# add_executable(geometry_dbg geometry_dbg.cpp)
# target_link_libraries(geometry_dbg)

# Debug grid
# add_executable(boundary_dbg boundary_dbg.cpp)
# target_link_libraries(boundary_dbg)
