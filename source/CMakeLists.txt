# CMakeList.txt : CMake project for main, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (3P22Lambda "main.cpp" "main.h")
target_include_directories(3P22Lambda PUBLIC ${PROJECT_BINARY_DIR})

# Add project version info
configure_file("config.h.in" "config.h" NEWLINE_STYLE LF)

# TODO: Add tests and install targets if needed.
