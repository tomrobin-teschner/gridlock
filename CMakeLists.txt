cmake_minimum_required(VERSION 3.23)

project(
  gridlock
  VERSION 0.9.0
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# List required dependencies here
find_package(Eigen3 REQUIRED)

# Add an executable and provide the source file directly.
add_executable(${CMAKE_PROJECT_NAME} "")

# Add nlohmann_json as a dependency
find_package(nlohmann_json REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} nlohmann_json::nlohmann_json)

# Add eigen as a dependency
find_package(Eigen3 REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} Eigen3::Eigen)

# Add ncurses as a dependency for fancy console printing
find_package(Curses REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} Curses::Curses)

# Traverse subdirectories and add source files synamically
add_subdirectory(src)

# Create a custom command to copy input files to the executable directory
add_custom_command(
    TARGET ${CMAKE_PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_SOURCE_DIR}/input" "${CMAKE_CURRENT_BINARY_DIR}/input"
)