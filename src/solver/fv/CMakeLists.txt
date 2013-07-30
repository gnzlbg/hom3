set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Unit Tests
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# # Heat Equation Solver:
# add_executable(heat_fv_solver_test heat_fv_solver_test.cpp)
# target_link_libraries(heat_fv_solver_test gtest gtest_main)
# add_test(heat_fv_solver_test heat_fv_solver_test)

# Euler Equation Solver:
# add_executable(euler_fv_solver_test euler_fv_solver_test.cpp)
# target_link_libraries(euler_fv_solver_test gtest gtest_main)
# add_test(euler_fv_solver_test euler_fv_solver_test)

# Coupled Heat-Euler Solver:
add_executable(coupled_heat_euler coupled_heat_euler_test.cpp)
target_link_libraries(coupled_heat_euler gtest gtest_main)
add_test(coupled_heat_euler coupled_heat_euler)

# Coupled Heat-Euler Solver:
add_executable(coupled_euler_euler coupled_euler_euler_test.cpp)
target_link_libraries(coupled_euler_euler gtest gtest_main)
add_test(coupled_euler_euler coupled_euler_euler)
