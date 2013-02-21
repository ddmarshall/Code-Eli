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
add_executable(FitContainerTest ${PROJECT_SOURCE_DIR}/fit_container/fit_container_test.cpp)
target_link_libraries(FitContainerTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME FitContainerTest COMMAND $<TARGET_FILE:FitContainerTest> --text)

# add the bezier_curve test
add_executable(BezierCurveTest ${PROJECT_SOURCE_DIR}/bezier_curve/bezier_curve_test.cpp)
target_link_libraries(BezierCurveTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME BezierCurveTest COMMAND $<TARGET_FILE:BezierCurveTest> --text)

# add the piecewise_curve test
add_executable(PiecewiseCurveTest ${PROJECT_SOURCE_DIR}/piecewise_curve/piecewise_curve_test.cpp)
target_link_libraries(PiecewiseCurveTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME PiecewiseCurveTest COMMAND $<TARGET_FILE:PiecewiseCurveTest> --text)

# add the explicit_bezier_curve test
add_executable(ExplicitBezierCurveTest ${PROJECT_SOURCE_DIR}/explicit_bezier_curve/explicit_bezier_curve_test.cpp)
target_link_libraries(ExplicitBezierCurveTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ExplicitBezierCurveTest COMMAND $<TARGET_FILE:ExplicitBezierCurveTest> --text)

#### add the cst_airfoil_curve test
###add_executable(CSTAirfoilCurveTest ${PROJECT_SOURCE_DIR}/cst_airfoil_curve/cst_airfoil_curve_test.cpp)
###target_link_libraries(CSTAirfoilCurveTest ${_ELI_TEST_LIBRARIES_LIST})
###add_test(NAME CSTAirfoilCurveTest COMMAND $<TARGET_FILE:CSTAirfoilCurveTest> --text)

# add the bezier_surface test
add_executable(BezierSurfaceTest ${PROJECT_SOURCE_DIR}/bezier_surface/bezier_surface_test.cpp)
target_link_libraries(BezierSurfaceTest ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME BezierSurfaceTest COMMAND $<TARGET_FILE:BezierSurfaceTest> --text)

