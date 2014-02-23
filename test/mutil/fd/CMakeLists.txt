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

project(FD)

# add the finite difference test
set (FD_TEST_EXE FDTest)
set (FD_TEST_SOURCE ${PROJECT_SOURCE_DIR}/fd_test.cpp)
set (FD_TEST_HEADER ${PROJECT_SOURCE_DIR}/fd_test_suite.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/fd.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/fd/d1o1.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/fd/d1o2.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/fd/d1o3.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/fd/d2o1.hpp
                    ${ELI_SOURCE_DIR}/include/eli/mutil/fd/d2o2.hpp)

add_executable(${FD_TEST_EXE} ${FD_TEST_SOURCE} ${FD_TEST_HEADER})
target_link_libraries(${FD_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${FD_TEST_EXE} COMMAND $<TARGET_FILE:${FD_TEST_EXE}> --text)
