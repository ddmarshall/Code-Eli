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

project(AD)

# set what libraries need to be used
set(_AD_LIBRARIES_LIST ${CPPTEST_LIBRARIES})
if (QD_FOUND)
  SET (_AD_LIBRARIES_LIST ${_AD_LIBRARIES_LIST} ${QD_LIBRARIES})
endif()

# add the mutil::constants test
add_executable(ADTest ${PROJECT_SOURCE_DIR}/test/ad_test.cpp)
target_link_libraries(ADTest ${_AD_LIBRARIES_LIST})
add_test(NAME ADTest COMMAND $<TARGET_FILE:ADTest> --text)
