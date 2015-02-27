cmake_minimum_required(VERSION 2.8)
project(trento CXX)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" CXX11)
if(NOT CXX11)
  message(FATAL_ERROR "Your compiler does not support C++11")
endif()

set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_PATCH dev)
set(PROJECT_VERSION
  "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")

# default build type: Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# default install prefix: ~/.local
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local"
      CACHE PATH "Install path prefix, prepended onto install directories."
      FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

# warnings adapted from http://stackoverflow.com/a/9862800
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -Wall -Wextra -Werror -Wpedantic -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winline -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=5 -Wundef")

set(LIBRARY_NAME "lib${PROJECT_NAME}")

add_subdirectory(src)
