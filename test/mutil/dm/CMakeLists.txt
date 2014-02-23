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

project(DM)

# add the discrete math test
set (DM_TEST_EXE DMTest)
set (DM_TEST_SOURCE ${PROJECT_SOURCE_DIR}/dm_test.cpp)
set (DM_TEST_HEADER ${PROJECT_SOURCE_DIR}/combination_test_suite.hpp
                    ${PROJECT_SOURCE_DIR}/factorial_test_suite.hpp
                    ${PROJECT_SOURCE_DIR}/binomial_coefficient_test_suite.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/dm/factorial.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/dm/combination.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/dm/binomial_coefficient.hpp)

add_executable(${DM_TEST_EXE} ${DM_TEST_SOURCE} ${DM_TEST_HEADER})
target_link_libraries(${DM_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${DM_TEST_EXE} COMMAND $<TARGET_FILE:${DM_TEST_EXE}> --text)
