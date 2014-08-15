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

project(ELI_TEST)

# include common test header files
set(_ELI_TEST_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
include_directories(${_ELI_TEST_INCLUDE_DIR})

# set what libraries need to be used
set(_ELI_TEST_LIBRARIES_LIST ${CPPTEST_LIBRARIES})

add_subdirectory(constants)
add_subdirectory(mutil)
add_subdirectory(geom)
add_subdirectory(util)
