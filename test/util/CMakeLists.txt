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
add_executable(ToleranceTest ${PROJECT_SOURCE_DIR}/tolerance_test.cpp)
target_link_libraries(ToleranceTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ToleranceTest COMMAND $<TARGET_FILE:ToleranceTest> --text)

# add the floating point test
add_executable(FloatingPointTest ${PROJECT_SOURCE_DIR}/floating_point_test.cpp)
target_link_libraries(FloatingPointTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME FloatingPointTest COMMAND $<TARGET_FILE:FloatingPointTest> --text)

# add the traits test
add_executable(TraitsTest ${PROJECT_SOURCE_DIR}/traits_test.cpp)
target_link_libraries(TraitsTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME TraitsTest COMMAND $<TARGET_FILE:TraitsTest> --text)
