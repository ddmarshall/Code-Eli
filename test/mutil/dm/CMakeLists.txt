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
add_executable(DMTest ${PROJECT_SOURCE_DIR}/dm_test.cpp)
target_link_libraries(DMTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME DMTest COMMAND $<TARGET_FILE:DMTest> --text)
