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

project(QUAD)

# add the quadrature test
set (QUAD_TEST_EXE QuadTest)
set (QUAD_TEST_SOURCE ${PROJECT_SOURCE_DIR}/quad_test.cpp)
set (QUAD_TEST_HEADER ${PROJECT_SOURCE_DIR}/trapezoid_test_suite.hpp
                      ${PROJECT_SOURCE_DIR}/simpson_test_suite.hpp
                      ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                      ${ELI_SOURCE_DIR}/include/eli/mutil/quad/trapezoid.hpp
                      ${ELI_SOURCE_DIR}/include/eli/mutil/quad/simpson.hpp)

add_executable(${QUAD_TEST_EXE} ${QUAD_TEST_SOURCE} ${QUAD_TEST_HEADER})
target_link_libraries(${QUAD_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${QUAD_TEST_EXE} COMMAND $<TARGET_FILE:${QUAD_TEST_EXE}> --text)
