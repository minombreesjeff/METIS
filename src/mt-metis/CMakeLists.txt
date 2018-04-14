cmake_minimum_required(VERSION 2.6)
project(MTMETIS)

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

if (DEFINED DEBUG OR DEFINED ASSERT)
  if (DEFINED DEBUG)
    add_definitions(-DDEBUG=${DEBUG})
    message("Debugging enabled")
  endif()
  if (DEFINED ASSERT)
    add_definitions(-DUSE_ASSERTS=${ASSERT})
    message("Assertions enabled")
  else()
    add_definitions(-DNDEBUG=1)
    add_definitions(-DNDEBUG2=1)
  endif()
  set(CMAKE_BUILD_TYPE Debug)
  message("Building a debug release")
else()
  set(CMAKE_BUILD_TYPE Release)
  add_definitions(-DNDEBUG=1)
  add_definitions(-DNDEBUG2=1)
  message("Building an optimized release")
endif()

if (SHARED)
  set(MTMETIS_LIBRARY_TYPE SHARED)
else()
  set(MTMETIS_LIBRARY_TYPE STATIC)
endif()

if (DEFINED BIGEDGES)
  add_definitions(-DBOWSTRING_64BIT_EDGES=${BIGEDGES})
  add_definitions(-DMTMETIS_64BIT_EDGES=${BIGEDGES})
endif()

if (DEFINED BIGVERTICES)
  add_definitions(-DBOWSTRING_64BIT_VERTICES=${BIGVERTICES})
  add_definitions(-DMTMETIS_64BIT_VERTICES=${BIGVERTICES})
endif()

if (DEFINED BIGWEIGHTS)
  add_definitions(-DBOWSTRING_64BIT_WEIGHTS=${BIGWEIGHTS})
  add_definitions(-DMTMETIS_64BIT_WEIGHTS=${BIGWEIGHTS})
endif()

if (DEFINED BIGEDGES OR DEFINED BIGVERTICES OR DEFINED BIGWEIGHTS)
  add_definitions(-DIDXTYPEWIDTH=64)
endif()

if (NOT DEFINED DISABLE_ZLIB OR NOT DISABLE_ZLIB)
  set(REQ_LIBS ${REQ_LIBS} z)
else()
  add_definitions(-DDISABLE_ZLIB=${DISABLE_ZLIB})
endif()
set(REQ_LIBS ${REQ_LIBS} m)

if (DEFINED DEVEL)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wmissing-prototypes") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wmissing-declarations") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wstrict-aliasing") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-parameter") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pedantic") 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fstrict-aliasing") 
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99") 

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -no-ansi-alias")
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# make sure we get the POSIX things
add_definitions(-D_POSIX_C_SOURCE=1)

# handle graph utils wgts
add_definitions(-DBOWSTRING_INT_WEIGHTS=1)

include_directories("${CMAKE_SOURCE_DIR}/include")
include_directories("${DOMLIB_PATH}")
include_directories("${BOWSTRING_PATH}/include")
include_directories("${METIS_PATH}/include")
include_directories("${METIS_PATH}/libmetis")
include_directories("${METIS_PATH}/GKlib")
set(BOWSTRING_LIB_PATH "${CMAKE_SOURCE_DIR}/${BOWSTRING_PATH}/src")
set(METIS_LIB_PATH "${CMAKE_SOURCE_DIR}/${METIS_PATH}/libmetis")
add_subdirectory("include")
add_subdirectory("src")

if (DEFINED TESTS)
  enable_testing()
  add_subdirectory("apitest")
endif()



