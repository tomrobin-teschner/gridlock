cmake_minimum_required(VERSION 3.23)

project(
  gridlock
  VERSION 0.10.0
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# List required dependencies here
find_package(Eigen3 REQUIRED)

# Add an executable and provide the source file directly.
add_executable(${CMAKE_PROJECT_NAME} "")

# Add tomlcplusplus as a dependency
find_package(tomlplusplus REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} tomlplusplus::tomlplusplus)

# Add eigen as a dependency
find_package(Eigen3 REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} Eigen3::Eigen)

# Add ncurses as a dependency for fancy console printing
find_package(Curses REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} Curses::Curses)

# Traverse subdirectories and add source files synamically
add_subdirectory(src)