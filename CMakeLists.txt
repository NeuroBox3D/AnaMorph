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
set(cxx_flags "${warnings} ${exclude_warnings} ${additionalFlags} -std=${cpp_standard} -pedantic -O3 -pthread")

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



