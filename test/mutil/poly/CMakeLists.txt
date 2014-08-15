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

project(POLY_TEST)

# add the polynomial test
set (POLY_TEST_EXE PolyTest)
set (POLY_TEST_SOURCE ${PROJECT_SOURCE_DIR}/poly_test.cpp)
set (POLY_TEST_HEADER ${PROJECT_SOURCE_DIR}/poly_basic_test_suite.hpp
                      ${PROJECT_SOURCE_DIR}/poly_math_test_suite.hpp
                      ${PROJECT_SOURCE_DIR}/poly_root_test_suite.hpp
                      ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                      ${ELI_SOURCE_DIR}/include/eli/mutil/dm/combination.hpp
                      ${ELI_SOURCE_DIR}/include/eli/mutil/poly/polynomial.hpp
                      ${ELI_SOURCE_DIR}/include/eli/mutil/poly/root/descartes_rule.hpp
                      ${ELI_SOURCE_DIR}/include/eli/mutil/poly/root/sturm_count.hpp
                      ${ELI_SOURCE_DIR}/include/eli/mutil/poly/root/closed_form.hpp
                      ${ELI_SOURCE_DIR}/include/eli/mutil/poly/root/radius.hpp
                      ${ELI_SOURCE_DIR}/include/eli/mutil/poly/root/sign_changes.hpp)

add_executable(${POLY_TEST_EXE} ${POLY_TEST_SOURCE} ${POLY_TEST_HEADER})
target_link_libraries(${POLY_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${POLY_TEST_EXE} COMMAND $<TARGET_FILE:${POLY_TEST_EXE}> --text)
