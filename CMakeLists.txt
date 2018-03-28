project(MPC)

cmake_minimum_required(VERSION 3.5)

add_definitions(-std=c++11 -o3)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(sources src/main.cpp src/mpc.cpp src/cost.cpp)

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

    include_directories(/usr/local/opt/openssl/include)
    link_directories(/usr/local/opt/openssl/lib)
    link_directories(/usr/local/Cellar/libuv/1*/lib)

endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

add_executable(mpc ${sources})

target_link_libraries(mpc ipopt z ssl uv uWS)

