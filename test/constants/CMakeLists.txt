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

project(CONSTANTS)

# add the constants test
set (CONSTANTS_MATH_TEST_EXE ConstantsMathTest)
set (CONSTANTS_MATH_TEST_SOURCE ${PROJECT_SOURCE_DIR}/constants_math_test.cpp)
set (CONSTANTS_MATH_TEST_HEADER ${PROJECT_SOURCE_DIR}/constants_math_test_suite.hpp
                                ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp)

add_executable(${CONSTANTS_MATH_TEST_EXE} ${CONSTANTS_MATH_TEST_SOURCE} ${CONSTANTS_MATH_TEST_HEADER})
target_link_libraries(${CONSTANTS_MATH_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${CONSTANTS_MATH_TEST_EXE} COMMAND $<TARGET_FILE:${CONSTANTS_MATH_TEST_EXE}> --text)
