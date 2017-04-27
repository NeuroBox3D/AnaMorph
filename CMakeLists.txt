##
 # This file is part of
 #
 # AnaMorph: a framework for geometric modelling, consistency analysis and surface
 # mesh generation of anatomically reconstructed neuron morphologies.
 # 
 # Copyright (c) 2013-2017: G-CSC, Goethe University Frankfurt - Queisser group
 # Author: Markus Breit
 # 
 # AnaMorph is free software: Redistribution and use in source and binary forms,
 # with or without modification, are permitted under the terms of the
 # GNU Lesser General Public License version 3 (as published by the
 # Free Software Foundation) with the following additional attribution
 # requirements (according to LGPL/GPL v3 §7):
 # 
 # (1) The following notice must be displayed in the Appropriate Legal Notices
 # of covered and combined works:
 # "Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph)."
 # 
 # (2) The following notice must be displayed at a prominent place in the
 # terminal output of covered works:
 # "Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph)."
 #
 # (3) Neither the name "AnaMorph" nor the names of its contributors may be
 # used to endorse or promote products derived from this software without
 # specific prior written permission.
 #
 # (4) The following bibliography is recommended for citation and must be
 # preserved in all covered files:
 # "Mörschel K, Breit M, Queisser G. Generating neuron geometries for detailed
 #   three-dimensional simulations using AnaMorph. Neuroinformatics (2017)"
 # "Grein S, Stepniewski M, Reiter S, Knodel MM, Queisser G.
 #   1D-3D hybrid modelling – from multi-compartment models to full resolution
 #   models in space and time. Frontiers in Neuroinformatics 8, 68 (2014)"
 # "Breit M, Stepniewski M, Grein S, Gottmann P, Reinhardt L, Queisser G.
 #   Anatomically detailed and large-scale simulations studying synapse loss
 #   and synchrony using NeuroBox. Frontiers in Neuroanatomy 10 (2016)"
 # 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 # See the GNU Lesser General Public License for more details.
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program. If not, see <http://www.gnu.org/licenses/>.
 ##
 
cmake_minimum_required (VERSION 2.7)
project (AnaMorph)

## options
option(DEBUG "compile in debug mode" OFF)
message(STATUS "DEBUG     ${DEBUG}")

option(CELLGEN "build target am_cellgen" ON)
message(STATUS "CELLGEN   ${CELLGEN}")

option(MESHSTAT "build target am_meshstat" ON)
message(STATUS "MESHSTAT  ${MESHSTAT}")

option(SHARED "build shared library" OFF)
message(STATUS "SHARED    ${SHARED}")


## check if boost is available
set(boost_cmp_flag)
find_package(Boost)
if (Boost_FOUND)
	# we only need headers
	include_directories(${Boost_INCLUDE_DIRS})
	
	# set compile flag
	set(boost_cmp_flag -DWITH_BOOST)
endif (Boost_FOUND)

## directories
get_filename_component(ROOT_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ROOT_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT_DIR}/lib)

## compile flags
set(cpp_standard c++11)
#set(c_standard c99)
set(stdlib libc++)

# clang does not like -pthread in the linking stage
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set(additionalFlags "${additional_flags} -Qunused-arguments")
endif()
set(exclude_warnings "-Wno-variadic-macros -Wno-missing-braces -Wno-format -Wno-unused-parameter")
set(warnings "-Wall -Wextra")
set(cxx_flags "${warnings} ${exclude_warnings} ${additionalFlags} ${boost_cmp_flag} -std=${cpp_standard} -pedantic -O3 -pthread")

if (DEBUG)
	set(cxx_flags "${cxx_flags} -g -D__DEBUG__")
else (DEBUG)
	set(cxx_flags "${cxx_flags} -DNDEBUG")
endif (DEBUG)

set(CMAKE_CXX_FLAGS "" CACHE STRING "clear flags" FORCE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${cxx_flags}" CACHE STRING "overriden flags!" FORCE)


## include dirs
include_directories(include tsrc)

## build the library
if (SHARED)
	set(BUILD_SHARED_LIBS ON)
else (SHARED)
	set(BUILD_SHARED_LIBS OFF)
endif (SHARED)
if (${CMAKE_VERSION} VERSION_GREATER 3.0
	OR ${CMAKE_VERSION} VERSION_EQUAL 3.0)
	cmake_policy(SET CMP0042 NEW)
endif()

set(AMLIB_SOURCES
	src/aux.cc
	src/IdQueue.cc
	src/CLApplication.cc
	src/AnaMorph_cellgen.cc
	src/Vec3.cc
	src/Vec2.cc
)
add_library(anamorph ${AMLIB_SOURCES})

## targets
if (CELLGEN)	
	add_executable(am_cellgen src/am_cellgen.cc)
	target_link_libraries(am_cellgen anamorph) 
endif (CELLGEN)

if (MESHSTAT)
	add_executable(am_meshstat src/am_meshstat.cc)
	target_link_libraries(am_meshstat anamorph)
endif (MESHSTAT)



