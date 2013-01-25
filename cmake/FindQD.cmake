# Find QD library
#
# This tries to find the QD library. Once completed it will set:
# QD_INCLUDE_DIRS - include directories for compilation
# QD_LIBRARIES    - libraries needed to link against
# QD_LIBRARY_DIRS - directories where to find the libraries
# QD_FOUND        - boolean set to true if found
#
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

if(NOT QD_FIND_VERSION)
  if(NOT QD_FIND_VERSION_MAJOR)
    set(QD_FIND_VERSION_MAJOR 0)
  endif()
  if(NOT QD_FIND_VERSION_MINOR)
    set(QD_FIND_VERSION_MINOR 0)
  endif()
  if(NOT QD_FIND_VERSION_PATCH)
    set(QD_FIND_VERSION_PATCH 0)
  endif()

  set(QD_FIND_VERSION "${QD_FIND_VERSION_MAJOR}.${QD_FIND_VERSION_MINOR}.${QD_FIND_VERSION_PATCH}")
endif()

macro(_qd_check_version)
  if ((QD_VERSION VERSION_GREATER QD_FIND_VERSION) OR (QD_VERSION VERSION_EQUAL QD_FIND_VERSION) )
    set(QD_VERSION_OK TRUE)
  else()
    set(QD_VERSION_OK FALSE)
    STRING(REGEX REPLACE "(\r?\n)+$" "" QD_VERSION_STR "${QD_VERSION}")
    message(STATUS "Found QD version ${QD_VERSION_STR}, but version ${QD_FIND_VERSION} is required")
  endif()
endmacro()

if(QD_LIBRARIES AND QD_LIBRARIES AND QD_LIBRARY_DIRS)
  _qd_check_version()
else()

  find_program(QD_CONFIG qd-config
               PATHS
              /usr/bin
              /usr/local/bin
              /opt/local/bin
              ${CMAKE_INSTALL_PREFIX}/bin
              $ENV{QDDIR}
              $ENV{QDDIR}/bin)

  if (NOT (QD_CONFIG STREQUAL "QD_CONFIG-NOTFOUND"))
    execute_process(COMMAND ${QD_CONFIG} --version OUTPUT_VARIABLE QD_VERSION)

    _qd_check_version()

    if (QD_VERSION_OK)
      find_path(QD_INCLUDE_DIRS qd/qd_real.h
                PATHS
                  /usr/include
                  /usr/local/include
                  /opt/local/include
                  ${CMAKE_INSTALL_PREFIX}/include
                  $ENV{QDDIR}
                  ${INCLUDE_INSTALL_DIR})
      find_library(QD_LIBRARIES qd
                   PATHS
                     $ENV{QDDIR}
                     /usr/local/lib
                     /opt/local/lib
                     ${LIB_INSTALL_DIR})
      get_filename_component(QD_LIBRARY_DIRS ${QD_LIBRARIES} PATH)
    endif()
  endif()

  mark_as_advanced(QD_INCLUDE_DIRS QD_LIBRARIES QD_FOUND)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QD DEFAULT_MSG QD_INCLUDE_DIRS QD_LIBRARIES QD_LIBRARY_DIRS QD_VERSION_OK)
