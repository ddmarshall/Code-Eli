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

# add the bezier_surface test
set (BEZIER_SURFACE_TEST_EXE BezierSurfaceTest)
set (BEZIER_SURFACE_TEST_SOURCE ${PROJECT_SOURCE_DIR}/bezier_surface_test.cpp)
set (BEZIER_SURFACE_TEST_HEADER ${PROJECT_SOURCE_DIR}/bezier_surface_test_suite.hpp
                                ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                                ${ELI_SOURCE_DIR}/include/eli/util/tolerance.hpp
                                ${ELI_SOURCE_DIR}/include/eli/geom/general/continuity.hpp
                                ${ELI_SOURCE_DIR}/include/eli/geom/general/bounding_box.hpp
                                ${ELI_SOURCE_DIR}/include/eli/geom/point/distance.hpp
                                ${ELI_SOURCE_DIR}/include/eli/geom/curve/bezier.hpp
                                ${ELI_SOURCE_DIR}/include/eli/geom/curve/equivalent_curves.hpp
                                ${ELI_SOURCE_DIR}/include/eli/geom/surface/bezier.hpp
                                ${ELI_SOURCE_DIR}/include/eli/geom/surface/curvature.hpp)

add_executable(${BEZIER_SURFACE_TEST_EXE} ${BEZIER_SURFACE_TEST_SOURCE} ${BEZIER_SURFACE_TEST_HEADER})
target_link_libraries(${BEZIER_SURFACE_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${BEZIER_SURFACE_TEST_EXE} COMMAND $<TARGET_FILE:${BEZIER_SURFACE_TEST_EXE}> --text)
