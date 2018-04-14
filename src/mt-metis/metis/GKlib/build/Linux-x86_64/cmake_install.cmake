# Install script for directory: /home/grad03/lasalle/src/GKlib

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/grad03/lasalle/src/GKlib/build/Linux-x86_64/libGKlib.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/grad03/lasalle/src/GKlib/gk_arch.h"
    "/home/grad03/lasalle/src/GKlib/gk_macros.h"
    "/home/grad03/lasalle/src/GKlib/GKlib.h"
    "/home/grad03/lasalle/src/GKlib/gk_struct.h"
    "/home/grad03/lasalle/src/GKlib/gk_types.h"
    "/home/grad03/lasalle/src/GKlib/gk_mkpqueue.h"
    "/home/grad03/lasalle/src/GKlib/gk_mkmemory.h"
    "/home/grad03/lasalle/src/GKlib/gk_defs.h"
    "/home/grad03/lasalle/src/GKlib/ms_stat.h"
    "/home/grad03/lasalle/src/GKlib/ms_stdint.h"
    "/home/grad03/lasalle/src/GKlib/gkregex.h"
    "/home/grad03/lasalle/src/GKlib/gk_mkrandom.h"
    "/home/grad03/lasalle/src/GKlib/gk_mkblas.h"
    "/home/grad03/lasalle/src/GKlib/gk_mksort.h"
    "/home/grad03/lasalle/src/GKlib/gk_proto.h"
    "/home/grad03/lasalle/src/GKlib/gk_mkutils.h"
    "/home/grad03/lasalle/src/GKlib/gk_externs.h"
    "/home/grad03/lasalle/src/GKlib/ms_inttypes.h"
    "/home/grad03/lasalle/src/GKlib/gk_getopt.h"
    "/home/grad03/lasalle/src/GKlib/gk_mkpqueue2.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/grad03/lasalle/src/GKlib/build/Linux-x86_64/test/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/home/grad03/lasalle/src/GKlib/build/Linux-x86_64/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/home/grad03/lasalle/src/GKlib/build/Linux-x86_64/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
