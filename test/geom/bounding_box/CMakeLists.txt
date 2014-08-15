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

project(BOUNDING_BOX_TEST)

#add the bounding box test
set (BOUNDING_BOX_TEST_EXE BoundingBoxTest)
set (BOUNDING_BOX_TEST_SOURCE ${PROJECT_SOURCE_DIR}/bounding_box_test.cpp)
set (BOUNDING_BOX_TEST_HEADER ${PROJECT_SOURCE_DIR}/bounding_box_test_suite.hpp
                              ${ELI_SOURCE_DIR}/include/eli/util/tolerance.hpp
                              ${ELI_SOURCE_DIR}/include/eli/geom/general/bounding_box.hpp)

add_executable(${BOUNDING_BOX_TEST_EXE} ${BOUNDING_BOX_TEST_SOURCE} ${BOUNDING_BOX_TEST_HEADER})
target_link_libraries(${BOUNDING_BOX_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${BOUNDING_BOX_TEST_EXE} COMMAND $<TARGET_FILE:${BOUNDING_BOX_TEST_EXE}> --text)
