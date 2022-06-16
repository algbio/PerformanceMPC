add_executable(wac wac.cpp)
target_link_libraries(wac mpc_lib nlohmann_json::nlohmann_json CLI11::CLI11)

