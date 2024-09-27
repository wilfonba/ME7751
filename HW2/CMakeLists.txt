# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(main VERSION 1.0)

# Specify the C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Automatically find all .c files in the src/ directory
file(GLOB SOURCES
    ${PROJECT_SOURCE_DIR}/src/*.c
)

# Find the libxml2 package
find_package(LibXml2 REQUIRED)

# Create the executable
add_executable(main ${SOURCES})

# Link with libxml2
target_link_libraries(main ${LIBXML2_LIBRARIES})

# Include the libxml2 headers
target_include_directories(main PRIVATE ${LIBXML2_INCLUDE_DIR})
