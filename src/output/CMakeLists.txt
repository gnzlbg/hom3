set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Unit Tests
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# File Handler Test:
#add_executable(file_test file_test.cpp)
# target_link_libraries(file_test gtest gtest_main ${Boost_LIBRARIES}
#                       ${MPI_LIBRARIES} ${PNETCDF_LIBRARY})
#target_link_libraries(file_test gtest gtest_main ${HOM3_BASE_LIBS})
#add_test(file_test file_test)
