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

project(NLS)

# set what libraries need to be used
set(_NLS_LIBRARIES_LIST ${CPPTEST_LIBRARIES})
if (QD_FOUND)
  SET (_NLS_LIBRARIES_LIST ${_NLS_LIBRARIES_LIST} ${QD_LIBRARIES})
endif()

# add the nonlinear solver test
add_executable(NLSTest ${PROJECT_SOURCE_DIR}/test/nls_test.cpp)
target_link_libraries(NLSTest ${_NLS_LIBRARIES_LIST})
add_test(NAME NLSTest COMMAND $<TARGET_FILE:NLSTest> --text)
