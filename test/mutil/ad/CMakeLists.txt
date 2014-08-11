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

project(AD_TEST)

# add the automatic differentiation test
set (AD_TEST_EXE ADTest)
set (AD_TEST_SOURCE ${PROJECT_SOURCE_DIR}/ad_test.cpp)
set (AD_TEST_HEADER ${PROJECT_SOURCE_DIR}/ad_dual_test_suite.hpp
                    ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                    ${ELI_SOURCE_DIR}/include/eli/util/traits.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/ad/dual.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/ad/dual_number.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/ad/dual_functions.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/ad/dual_operators.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/ad/dual_math.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/ad/dual_traits.hpp)

add_executable(${AD_TEST_EXE} ${AD_TEST_SOURCE} ${AD_TEST_HEADER})
target_link_libraries(${AD_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${AD_TEST_EXE} COMMAND $<TARGET_FILE:${AD_TEST_EXE}> --text)
