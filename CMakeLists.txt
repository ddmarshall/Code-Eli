################################################################################
# Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
#
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the Eclipse Public License v1.0
# which accompanies this distribution, and is available at
# http://www.eclipse.org/legal/epl-v10.html
#
# Contributors:
#    David D. Marshall - initial code and implementation
################################################################################

cmake_minimum_required(VERSION 2.8)

project(ELI)

# set the release version number
set(ELI_VERSION_MAJOR 0)
set(ELI_VERSION_MINOR 3)
set(ELI_VERSION_PATCH 0)

# capture the build date and time
if (WIN32)
  execute_process(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE _ELI_WIN_DATE)
  string(REGEX MATCH "(..)/(..)/(....)" _ELI_WIN_DATE ${_ELI_WIN_DATE})
  string(REGEX REPLACE "(..)/(..)/(....)" "\\3\\1\\2" ELI_BUILD_DATE ${_ELI_WIN_DATE})
  execute_process(COMMAND "cmd" " /C echo %time%" OUTPUT_VARIABLE _ELI_WIN_TIME)
  string(REGEX MATCH "(..):(..):(..)" _ELI_WIN_TIME ${_ELI_WIN_TIME})
  string(REGEX REPLACE "(..):(..):(..)" "\\1" _ELI_WIN_HOUR_TIME ${_ELI_WIN_TIME})
  string(STRIP ${_ELI_WIN_HOUR_TIME} _ELI_WIN_HOUR_TIME)
  if(_ELI_WIN_HOUR_TIME LESS "10")
    set(_ELI_WIN_HOUR_TIME "0${_ELI_WIN_HOUR_TIME}")
  endif()
  string(REGEX REPLACE "(..):(..):(..)" "\\2\\3" _ELI_WIN_MINSEC_TIME ${_ELI_WIN_TIME})
  set(ELI_BUILD_TIME "${_ELI_WIN_HOUR_TIME}${_ELI_WIN_MINSEC_TIME}")
elseif (UNIX)
  execute_process(COMMAND "date" "+%Y%m%d" OUTPUT_VARIABLE ELI_BUILD_DATE)
  string(REGEX REPLACE "(\r?\n)+$" "" ELI_BUILD_DATE "${ELI_BUILD_DATE}")
  execute_process(COMMAND "date" "+%H%M%S" OUTPUT_VARIABLE ELI_BUILD_TIME)
  string(REGEX REPLACE "(\r?\n)+$" "" ELI_BUILD_TIME "${ELI_BUILD_TIME}")
else()
  message(SEND_ERROR "Could not determine build date and time.")
endif()

# add the Code-Eli specific cmake files to module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${ELI_SOURCE_DIR}/cmake)

# include the compiler specific settings
include(${ELI_SOURCE_DIR}/cmake/ConfigureCompiler.cmake)

# include the CTest framework
include(CTest)

# find libraries that are needed or optional
# find the cpptest library
find_package(Doxygen)
if (DOXYGEN_FOUND)
  configure_file(${ELI_SOURCE_DIR}/cmake/Doxyfile.in ${ELI_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
                    COMMAND ${DOXYGEN_EXECUTABLE} ${ELI_BINARY_DIR}/Doxyfile
                    SOURCES ${ELI_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${ELI_BINARY_DIR}
                    COMMENT "Generating API documentation with Doxygen")
endif()

# configure the main eli.hpp file now that know what packages and features
# are to be included
configure_file(${ELI_SOURCE_DIR}/cmake/code_eli.hpp.in ${ELI_BINARY_DIR}/include/eli/code_eli.hpp @ONLY)

find_package(CPPTest)
find_package(Eigen3 3.0.0)

if (EIGEN3_FOUND)
  # add the include directories and libraries that are needed
  include_directories(${ELI_BINARY_DIR}/include
                      ${ELI_SOURCE_DIR}/include
                      ${EIGEN3_INCLUDE_DIR})
  add_subdirectory(example)

  if(CPPTEST_FOUND)
    # add the CPPTest include and link directories
    include_directories(${CPPTEST_INCLUDE_DIRS})
    link_directories(${CPPTEST_LIBRARY_DIRS})
    add_subdirectory(test)
  endif()
endif()
