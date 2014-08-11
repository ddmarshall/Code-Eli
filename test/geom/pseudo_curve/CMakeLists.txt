################################################################################
# Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
#
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the Eclipse Public License v1.0
# which accompanies this distribution, and is available at
# http://www.eclipse.org/legal/epl-v10.html
#
# Contributors:
#    David D. Marshall - initial code and implementation
################################################################################

project(AIFOIL_TEST)

# add the airfoil test
set (AIRFOIL_TEST_EXE AirfoilTest)
set (AIRFOIL_TEST_SOURCE ${PROJECT_SOURCE_DIR}/airfoil_test.cpp)
set (AIRFOIL_TEST_HEADER ${PROJECT_SOURCE_DIR}/four_digit_test_suite.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/util/tolerance.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/mutil/fd/d1o2.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/mutil/fd/d2o2.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/curve/pseudo/four_digit.hpp)

add_executable(${AIRFOIL_TEST_EXE} ${AIRFOIL_TEST_SOURCE} ${AIRFOIL_TEST_HEADER})
target_link_libraries(${AIRFOIL_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${AIRFOIL_TEST_EXE} COMMAND $<TARGET_FILE:${AIRFOIL_TEST_EXE}> --text)
