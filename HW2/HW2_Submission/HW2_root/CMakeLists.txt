# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(MyFortranProject LANGUAGES Fortran)

# Set the Fortran standard (optional)
set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Specify the source directory
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

# Find all Fortran source files in the src directory
file(GLOB_RECURSE Fortran_SOURCES "${SOURCE_DIR}/*.f90" "${SOURCE_DIR}/*.f" "${SOURCE_DIR}/*.f95")

# Add the executable target
add_executable(main ${Fortran_SOURCES})

# Link any additional libraries if needed (optional)
# target_link_libraries(MyFortranExecutable PRIVATE some_library)

# Specify include directories if needed (optional)
# target_include_directories(MyFortranExecutable PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Set compiler flags if needed (optional)
 set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")

# Add additional configurations as needed

# Install the executable (optional)
# install(TARGETS MyFortranExecutable DESTINATION bin)

