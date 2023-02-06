# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.24)

project (3P22Lambda
	VERSION 0.0.0.0
	DESCRIPTION "CRT phosphor decay simulator"
	HOMEPAGE_URL "https://github.com/Gumball2415/3P22Lambda"
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

# add tclap dependency
include_directories (BEFORE "source/tclap/include")

# Add project version info
configure_file("config.h.in" "config.h" NEWLINE_STYLE LF)

add_executable (3P22Lambda "source/main.cpp" "source/main.h")

# Add lodepng dependency
add_subdirectory ("source/lodepng")
list (APPEND EXTRA_LIBS "lodepng")
list (APPEND EXTRA_INCLUDES "${PROJECT_BINARY_DIR}/lodepng")

target_link_libraries (3P22Lambda PUBLIC ${EXTRA_LIBS})

target_include_directories (3P22Lambda
	PUBLIC "${PROJECT_BINARY_DIR}"
	${EXTRA_INCLUDES})

# TODO: Add tests and install targets if needed.