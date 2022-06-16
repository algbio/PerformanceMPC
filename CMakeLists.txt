cmake_minimum_required(VERSION 3.14)

project(PracticalMPC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE2 "-O3 -march=native -g -DNDEBUG")
#-DCMAKE_BUILD_TYPE=Release
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g ")
set(CMAKE_CXX_FLAGS_DEBUG2 "-O0 -g -fsanitize=address")

# Lemon and Eigen conflict
set_property(GLOBAL PROPERTY ALLOW_DUPLICATE_CUSTOM_TARGETS 1)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/15460959cbbfa20e66ef0b5ab497367e47fc0a04.zip
)
FetchContent_MakeAvailable(googletest)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.10.5/json.tar.xz)
FetchContent_MakeAvailable(json)

FetchContent_Declare(
  cli11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11
  GIT_TAG        291c58789c031208f08f4f261a858b5b7083e8e2
)
FetchContent_MakeAvailable(cli11)

FetchContent_Declare(
      Eigen
      GIT_REPOSITORY    https://gitlab.com/libeigen/eigen.git
      GIT_TAG           3.4.0
    )
FetchContent_MakeAvailable(Eigen)

include(FetchContent)
FetchContent_Declare(
  lemon
  URL file://${CMAKE_SOURCE_DIR}/ext/lemon-1.3.1.tar.gz
)
FetchContent_MakeAvailable(lemon)

include_directories(src)
include_directories(${LEMON_SOURCE_DIR})
include_directories(${LEMON_BINARY_DIR})


add_subdirectory("src")
add_subdirectory("test")
add_subdirectory("prog")
