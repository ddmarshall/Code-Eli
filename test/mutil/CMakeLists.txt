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

project(MUTIL)

if ( (CMAKE_CXX_COMPILER_ID STREQUAL "Intel") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "13.0") )
    message(WARNING "Version ${CMAKE_CXX_COMPILER_VERSION} of Intel compiler does not support the AD component.")
else()
  add_subdirectory(ad)
endif()
add_subdirectory(dm)
add_subdirectory(fd)
add_subdirectory(nls)
add_subdirectory(opt)
add_subdirectory(poly)
add_subdirectory(quad)
add_subdirectory(tolerance)
