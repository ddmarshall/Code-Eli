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

project(GEOM)

# add the fit_container test
set (FIT_CONTAINER_TEST_EXE FitContainerTest)
set (FIT_CONTAINER_TEST_SOURCE ${PROJECT_SOURCE_DIR}/fit_container_test.cpp)
set (FIT_CONTAINER_TEST_HEADER ${PROJECT_SOURCE_DIR}/fit_container_test_suite.hpp
                               ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                               ${ELI_SOURCE_DIR}/include/eli/mutil/fd/d1o2.hpp
                               ${ELI_SOURCE_DIR}/include/eli/mutil/fd/d2o2.hpp
                               ${ELI_SOURCE_DIR}/include/eli/geom/general/continuity.hpp
                               ${ELI_SOURCE_DIR}/include/eli/geom/point/distance.hpp
                               ${ELI_SOURCE_DIR}/include/eli/geom/curve/fit_container.hpp)

add_executable(${FIT_CONTAINER_TEST_EXE} ${FIT_CONTAINER_TEST_SOURCE} ${FIT_CONTAINER_TEST_HEADER})
target_link_libraries(${FIT_CONTAINER_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${FIT_CONTAINER_TEST_EXE} COMMAND $<TARGET_FILE:${FIT_CONTAINER_TEST_EXE}> --text)
