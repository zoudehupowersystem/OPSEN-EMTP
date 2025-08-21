cmake_minimum_required(VERSION 3.16)
project(EMTP_Simple_Case CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找Eigen库
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 添加头文件目录
include_directories(src)

# 添加可执行文件及其源文件
add_executable(emtp_sim
    main.cpp
    Control.cpp
    Grid.cpp
    Simulation.cpp
    Curve.cpp
    Devices.cpp
)

# 链接Eigen库
target_link_libraries(emtp_sim PRIVATE Eigen3::Eigen)

# 设置输出目录 (可选)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)