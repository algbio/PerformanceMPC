enable_testing()

add_executable(test test.cpp)

target_link_libraries(test PRIVATE mpc_lib)
target_link_libraries(test PRIVATE gtest_main)
target_link_libraries(test PRIVATE lemon)

include(GoogleTest)
gtest_discover_tests(test WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test)
