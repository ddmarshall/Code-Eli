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

project(POLY)

# set what libraries need to be used
set(_POLY_LIBRARIES_LIST ${CPPTEST_LIBRARIES})
if (QD_FOUND)
  SET (_POLY_LIBRARIES_LIST ${_POLY_LIBRARIES_LIST} ${QD_LIBRARIES})
endif()

# add the polynomial test
add_executable(PolyTest ${PROJECT_SOURCE_DIR}/test/poly_test.cpp)
target_link_libraries(PolyTest ${_POLY_LIBRARIES_LIST})
add_test(NAME PolyTest COMMAND $<TARGET_FILE:PolyTest> --text)
