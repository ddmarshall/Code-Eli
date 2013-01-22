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
set(ELI_VERSION_MINOR 0)
set(ELI_VERSION_PATCH 1)

# capture the build date and time
if (WIN32)
  execute_process(COMMAND "date" "/T" OUTPUT_VARIABLE _ELI_WIN_DATE)
  string(REGEX REPLACE "(..)/(..)/(....).*" "\\3\\1\\2" ELI_BUILD_DATE ${_ELI_WIN_DATE})
  execute_process(COMMAND "time" "/T" OUTPUT_VARIABLE _ELI_WIN_TIME)
  string(REGEX REPLACE "(.):(..):(..).*(.)" "\\1" _ELI_WIN_HOUR_TIME ${_ELI_WIN_TIME})
  string(REGEX REPLACE "(.):(..):(..).*(.)" "\\4" _ELI_WIN_AMPM_TIME ${_ELI_WIN_TIME})
  if(_ELI_WIN_AMPM_TIME STREQUAL "p")
    MATH(EXPR _ELI_WIN_HOUR_TIME "${_ELI_WIN_HOUR_TIME} + 12")
  else()
    if(_ELI_WIN_HOUR_TIME LESS "10")
      set(_ELI_WIN_HOUR_TIME "0${_ELI_WIN_HOUR_TIME}")
    endif()
  endif()
  string(REGEX REPLACE "(.):(..):(..).*(.)" "\\2\\3" _ELI_WIN_MINSEC_TIME ${_ELI_WIN_TIME})
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
find_package(CPPTest REQUIRED)
find_package(Eigen3 3.0.0 REQUIRED)
find_package(QD)
if (QD_FOUND)
  set(ELI_QD_FOUND 1)
endif ()

# configure the main eli.hpp file now that know what packages and features
# are to be included
configure_file(${ELI_SOURCE_DIR}/cmake/code_eli.hpp.in ${ELI_BINARY_DIR}/include/code_eli.hpp)

# add the include and link directories and libraries that are needed
include_directories(${ELI_BINARY_DIR}/include
                    ${PROJECT_SOURCE_DIR}/include
                    ${CPPTEST_INCLUDE_DIRS}
                    ${EIGEN3_INCLUDE_DIR})
link_directories(${CPPTEST_LIBRARY_DIRS})
if (QD_FOUND)
  include_directories(${QD_INCLUDE_DIRS})
  link_directories(${QD_LIBRARY_DIRS})
endif()

# add subdirectories of project
ADD_SUBDIRECTORY(constants)

