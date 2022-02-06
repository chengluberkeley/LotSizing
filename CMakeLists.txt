cmake_minimum_required (VERSION 3.0)
project(LotSizing)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
set(CMAKE_CXX_FLAGS "-std=c++11 -O3")

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/src
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

file(GLOB srcs
  ${PROJECT_SOURCE_DIR}/src/*.cpp
)

add_library(dynamic_path STATIC IMPORTED)
set_target_properties(dynamic_path PROPERTIES
  IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/DynamicPath/lib/libdynamic_path.a"
  INTERFACE_INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/DynamicPath/src"
)

add_executable(test_main ${PROJECT_SOURCE_DIR}/main.cpp ${srcs})
target_link_libraries(test_main PRIVATE dynamic_path)
