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

project(OPT_TEST)

# add the optimization test
set (OPT_TEST_EXE OptTest)
set (OPT_TEST_SOURCE ${PROJECT_SOURCE_DIR}/opt_test.cpp)
set (OPT_TEST_HEADER ${PROJECT_SOURCE_DIR}/least_squares_test_suite.hpp
                     ${ELI_SOURCE_DIR}/include/eli/mutil/opt/least_squares.hpp)

add_executable(${OPT_TEST_EXE} ${OPT_TEST_SOURCE} ${OPT_TEST_HEADER})
target_link_libraries(${OPT_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${OPT_TEST_EXE} COMMAND $<TARGET_FILE:${OPT_TEST_EXE}> --text)
