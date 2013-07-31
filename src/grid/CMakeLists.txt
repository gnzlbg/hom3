set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Unit Tests
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Debug grid
#add_executable(grid_dbg grid_dbg.cpp)
#target_link_libraries(grid_dbg)
