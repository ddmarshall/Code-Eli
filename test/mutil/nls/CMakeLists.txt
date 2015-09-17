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

project(NLS_TEST)

# add the nonlinear solver test
set (NLS_TEST_EXE NLSTest)
set (NLS_TEST_SOURCE ${PROJECT_SOURCE_DIR}/nls_test.cpp)
set (NLS_TEST_HEADER ${PROJECT_SOURCE_DIR}/nls_test_suite.hpp
                     ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                     ${ELI_SOURCE_DIR}/include/eli/mutil/nls.hpp
                     ${ELI_SOURCE_DIR}/include/eli/mutil/nls/iterative_root_base.hpp
                     ${ELI_SOURCE_DIR}/include/eli/mutil/nls/iterative_root_base_constrained.hpp
                     ${ELI_SOURCE_DIR}/include/eli/mutil/nls/bisection_method.hpp
                     ${ELI_SOURCE_DIR}/include/eli/mutil/nls/newton_raphson_method.hpp
                     ${ELI_SOURCE_DIR}/include/eli/mutil/nls/secant_method.hpp
                     ${ELI_SOURCE_DIR}/include/eli/mutil/nls/iterative_system_root_base.hpp
                     ${ELI_SOURCE_DIR}/include/eli/mutil/nls/iterative_system_root_base_constrained.hpp
                     ${ELI_SOURCE_DIR}/include/eli/mutil/nls/newton_raphson_system_method.hpp)

add_executable(${NLS_TEST_EXE} ${NLS_TEST_SOURCE} ${NLS_TEST_HEADER})
target_link_libraries(${NLS_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${NLS_TEST_EXE} COMMAND $<TARGET_FILE:${NLS_TEST_EXE}> --text)
