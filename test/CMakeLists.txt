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

# set what libraries need to be used
set(_ELI_TEST_LIBRARIES_LIST ${CPPTEST_LIBRARIES})
if (QD_FOUND)
  SET (_ELI_TEST_LIBRARIES_LIST ${_ELI_TEST_LIBRARIES_LIST} ${QD_LIBRARIES})
endif()

add_subdirectory(ad)
add_subdirectory(constants)
add_subdirectory(dm)
add_subdirectory(fd)
add_subdirectory(nls)
add_subdirectory(opt)
add_subdirectory(poly)
add_subdirectory(quad)
