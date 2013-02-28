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

project(POLY)

# add the polynomial test
add_executable(PolyTest ${PROJECT_SOURCE_DIR}/poly_test.cpp)
target_link_libraries(PolyTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME PolyTest COMMAND $<TARGET_FILE:PolyTest> --text)
