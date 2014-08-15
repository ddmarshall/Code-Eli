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

project(BEZIER_CURVE_TESTS)

# add the bezier_curve test
set (BEZIER_CURVE_TEST_EXE BezierCurveTest)
set (BEZIER_CURVE_TEST_SOURCE ${PROJECT_SOURCE_DIR}/bezier_curve_test.cpp)
set (BEZIER_CURVE_TEST_HEADER ${PROJECT_SOURCE_DIR}/bezier_curve_test_suite.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/util/tolerance.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/mutil/quad.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/mutil/opt/least_squares.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/mutil/dm/binomial_coefficient.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/utility/bezier.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/general/continuity.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/general/bounding_box.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/point/distance.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/curve/bezier.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/curve/length.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/curve/curvature.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/curve/fit_container.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/intersect/minimum_distance_curve.hpp)

add_executable(${BEZIER_CURVE_TEST_EXE} ${BEZIER_CURVE_TEST_SOURCE} ${BEZIER_CURVE_TEST_HEADER})
target_link_libraries(${BEZIER_CURVE_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${BEZIER_CURVE_TEST_EXE} COMMAND $<TARGET_FILE:${BEZIER_CURVE_TEST_EXE}> --text)

# add the bezier_curve_fit test
set (BEZIER_CURVE_FIT_TEST_EXE BezierCurveFitTest)
set (BEZIER_CURVE_FIT_TEST_SOURCE ${PROJECT_SOURCE_DIR}/bezier_curve_fit_test.cpp)
set (BEZIER_CURVE_FIT_TEST_HEADER ${PROJECT_SOURCE_DIR}/bezier_curve_fit_test_suite.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/util/tolerance.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/mutil/quad.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/mutil/opt/least_squares.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/mutil/dm/binomial_coefficient.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/utility/bezier.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/general/continuity.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/general/bounding_box.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/point/distance.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/curve/bezier.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/curve/length.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/curve/curvature.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/curve/fit_container.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/intersect/minimum_distance_curve.hpp)

add_executable(${BEZIER_CURVE_FIT_TEST_EXE} ${BEZIER_CURVE_FIT_TEST_SOURCE} ${BEZIER_CURVE_FIT_TEST_HEADER})
target_link_libraries(${BEZIER_CURVE_FIT_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${BEZIER_CURVE_FIT_TEST_EXE} COMMAND $<TARGET_FILE:${BEZIER_CURVE_FIT_TEST_EXE}> --text)
