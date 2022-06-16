add_subdirectory(wac)

add_executable(exp exp.cpp)
target_link_libraries(exp mpc_lib nlohmann_json::nlohmann_json CLI11::CLI11 Eigen3::Eigen)

