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

project(UTIL)

# add the tolerance test
set (TOLERANCE_TEST_EXE ToleranceTest)
set (TOLERANCE_TEST_SOURCE ${PROJECT_SOURCE_DIR}/tolerance_test.cpp)
set (TOLERANCE_TEST_HEADER ${PROJECT_SOURCE_DIR}/tolerance_test_suite.hpp
                           ${ELI_SOURCE_DIR}/include/eli/util/tolerance.hpp)

add_executable(${TOLERANCE_TEST_EXE} ${TOLERANCE_TEST_SOURCE} ${TOLERANCE_TEST_HEADER})
target_link_libraries(${TOLERANCE_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${TOLERANCE_TEST_EXE} COMMAND $<TARGET_FILE:${TOLERANCE_TEST_EXE}> --text)

# add the floating point test
set (FLOATING_POINT_TEST_EXE FloatingPointTest)
set (FLOATING_POINT_TEST_SOURCE ${PROJECT_SOURCE_DIR}/floating_point_test.cpp)
set (FLOATING_POINT_TEST_HEADER ${PROJECT_SOURCE_DIR}/floating_point_test_suite.hpp
                                ${ELI_SOURCE_DIR}/include/eli/util/floating_point.hpp
                                ${ELI_SOURCE_DIR}/include/eli/util/traits.hpp)

add_executable(${FLOATING_POINT_TEST_EXE} ${FLOATING_POINT_TEST_SOURCE} ${FLOATING_POINT_TEST_HEADER})
target_link_libraries(${FLOATING_POINT_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${FLOATING_POINT_TEST_EXE} COMMAND $<TARGET_FILE:${FLOATING_POINT_TEST_EXE}> --text)

# add the traits test
set (TRAITS_TEST_EXE TraitsTest)
set (TRAITS_TEST_SOURCE ${PROJECT_SOURCE_DIR}/traits_test.cpp)
set (TRAITS_TEST_HEADER ${PROJECT_SOURCE_DIR}/traits_test_suite.hpp
                        ${ELI_SOURCE_DIR}/include/eli/util/traits.hpp)

add_executable(${TRAITS_TEST_EXE} ${TRAITS_TEST_SOURCE} ${TRAITS_TEST_HEADER})
target_link_libraries(${TRAITS_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${TRAITS_TEST_EXE} COMMAND $<TARGET_FILE:${TRAITS_TEST_EXE}> --text)
