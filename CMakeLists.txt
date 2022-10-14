# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project (3P22Lambda
	VERSION 0.0.0.0
	DESCRIPTION "CRT phosphor decay emulator"
	HOMEPAGE_URL "https://github.com/Gumball2415/3P22Lambda")

# Add project version info
configure_file("config.h.in" "config.h" NEWLINE_STYLE LF)

# Include sub-projects.
add_subdirectory ("source")
