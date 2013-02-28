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
add_executable(QuadTest ${PROJECT_SOURCE_DIR}/quad_test.cpp)
target_link_libraries(QuadTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME QuadTest COMMAND $<TARGET_FILE:QuadTest> --text)
