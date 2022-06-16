add_library(mpc_lib STATIC graph.cpp naive.cpp lemon.cpp utils.cpp cc.cpp reach.cpp pflow.cpp antichain.cpp preprocessing.cpp transitive.cpp)
target_link_libraries(mpc_lib PRIVATE lemon nlohmann_json::nlohmann_json Eigen3::Eigen)
