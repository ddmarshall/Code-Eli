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

# add the minimum_distance test
set (MINIMUM_DISTANCE_TEST_EXE MinimumDistanceTest)
set (MINIMUM_DISTANCE_TEST_SOURCE ${PROJECT_SOURCE_DIR}/minimum_distance_test.cpp)
set (MINIMUM_DISTANCE_TEST_HEADER ${PROJECT_SOURCE_DIR}/minimum_distance_point_test_suite.hpp
                                  ${PROJECT_SOURCE_DIR}/minimum_distance_bounding_box_test_suite.hpp
                                  ${PROJECT_SOURCE_DIR}/minimum_distance_line_test_suite.hpp
                                  ${PROJECT_SOURCE_DIR}/minimum_distance_plane_test_suite.hpp
                                  ${PROJECT_SOURCE_DIR}/minimum_distance_curve_test_suite.hpp
                                  ${PROJECT_SOURCE_DIR}/minimum_distance_surface_test_suite.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/util/tolerance.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/mutil/nls/newton_raphson_constrained_method.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/general/bounding_box.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/point/distance.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/curve/bezier.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/surface/bezier.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/intersect/minimum_distance_point.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/intersect/minimum_distance_bounding_box.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/intersect/minimum_distance_line.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/intersect/minimum_distance_plane.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/intersect/minimum_distance_curve.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/intersect/minimum_distance_surface.hpp)


add_executable(${MINIMUM_DISTANCE_TEST_EXE} ${MINIMUM_DISTANCE_TEST_SOURCE} ${MINIMUM_DISTANCE_TEST_HEADER})
target_link_libraries(${MINIMUM_DISTANCE_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${MINIMUM_DISTANCE_TEST_EXE} COMMAND $<TARGET_FILE:${MINIMUM_DISTANCE_TEST_EXE}> --text)
