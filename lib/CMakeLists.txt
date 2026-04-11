cmake_minimum_required(VERSION 3.16)

project(terrain-descriptors LANGUAGES CXX)

# option to indicate we are building the library
option(TD_BUILD_LIB "Build terrain-descriptors as a library" ON)

# use C++17
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# compiler warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# pick sources only necessary for the library
set(TD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../appMetrics)

set(SOURCES
	${TD_DIR}/core.cpp
	${TD_DIR}/heightfield.cpp
	${TD_DIR}/heightfield-breach.cpp
	${TD_DIR}/heightfield-fill.cpp
	${TD_DIR}/heightfield-metrics-hydrology.cpp
	${TD_DIR}/heightfield-metrics-landforms.cpp
	${TD_DIR}/heightfield-metrics-local.cpp
	${TD_DIR}/heightfield-metrics-orometry.cpp
	${TD_DIR}/heightfield-metrics-visibility.cpp
	${TD_DIR}/integerfield.cpp
	${TD_DIR}/isolation.cpp
	${TD_DIR}/ppa.cpp
	${TD_DIR}/riversnet.cpp
	${TD_DIR}/scalarfield.cpp
	${TD_DIR}/terrainanalysis.cpp
	#
	${TD_DIR}/akirmse-mountains/divide_tree.cpp
	${TD_DIR}/akirmse-mountains/domain_map.cpp
	${TD_DIR}/akirmse-mountains/island_tree.cpp
	${TD_DIR}/akirmse-mountains/latlng.cpp
	${TD_DIR}/akirmse-mountains/line_tree.cpp
	${TD_DIR}/akirmse-mountains/tile.cpp
	${TD_DIR}/akirmse-mountains/tree_builder.cpp
	${TD_DIR}/akirmse-mountains/utm_coordinate_system.cpp
)

set(INCLUDE_DIR ${TD_DIR} ${TD_DIR}/akirmse-mountains)

add_library(${PROJECT_NAME} STATIC ${SOURCES})

# hacky way to force appMetrics as a prefix for the header includes
target_include_directories(
  ${PROJECT_NAME}
  PUBLIC
	$<BUILD_INTERFACE:${TD_DIR}/..>
  PRIVATE
	${TD_DIR}
	${TD_DIR}/akirmse-mountains	
)

# add preprocessor definition if building as a library
if(TD_BUILD_LIB)
    target_compile_definitions(${PROJECT_NAME} PUBLIC TD_BUILD_LIB)
endif()
