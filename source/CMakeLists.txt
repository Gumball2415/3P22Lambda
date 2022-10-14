# CMakeList.txt : CMake project for main, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (3P22Lambda "main.cpp" "main.h")

target_include_directories(3P22Lambda
	PUBLIC "${PROJECT_BINARY_DIR}"
	"${PROJECT_BINARY_DIR}/lodepng")

# Add lodepng dependency
add_subdirectory ("lodepng")
target_link_libraries (3P22Lambda PUBLIC lodepng)

# TODO: Add tests and install targets if needed.
