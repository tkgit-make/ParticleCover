cmake_minimum_required(VERSION 3.10)

project(ParticleCover VERSION 0.0.1)

# cpp standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Add subdirectories
add_subdirectory(src)

# Define the executable
add_executable(ParticleCover ${SOURCES})

# Include directories
target_include_directories(ParticleCover PUBLIC include)
