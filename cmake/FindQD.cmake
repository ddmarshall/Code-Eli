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

if(QD_LIBRARIES)
  set(QD_FOUND TRUE)
else()
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

  mark_as_advanced(QD_INCLUDE_DIRS QD_LIBRARIES QD_FOUND)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(QD DEFAULT_MSG QD_INCLUDE_DIRS QD_LIBRARIES QD_LIBRARY_DIRS)
endif()
