cmake_minimum_required(VERSION 3.24)

option(GALAXY_ENABLE_CUDA "Build CUDA parts of the engine" OFF)

if(GALAXY_ENABLE_CUDA)
    project(GalaxyEngine LANGUAGES CXX CUDA)
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
else()
    project(GalaxyEngine LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
if(NOT Python_EXECUTABLE AND DEFINED ENV{Python_EXECUTABLE})
    set(Python_EXECUTABLE $ENV{Python_EXECUTABLE})
endif()
message(STATUS "Building engine for Python ${Python_VERSION}")
message(STATUS "  executable : ${Python_EXECUTABLE}")
message(STATUS "  include    : ${Python_INCLUDE_DIRS}")

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

add_subdirectory(src/engine)
add_subdirectory(src/bindings)

set(ENGINE_OUTPUT_DIR "${CMAKE_SOURCE_DIR}" CACHE PATH "Where to copy the engine module")