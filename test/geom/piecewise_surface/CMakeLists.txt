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

# add the piecewise_surface test
set (PIECEWISE_SURFACE_TEST_EXE PiecewiseSurfaceTest)
set (PIECEWISE_SURFACE_TEST_SOURCE ${PROJECT_SOURCE_DIR}/piecewise_surface_test.cpp)
set (PIECEWISE_SURFACE_TEST_HEADER ${PROJECT_SOURCE_DIR}/piecewise_surface_test_suite.hpp
                                   ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                                   ${ELI_SOURCE_DIR}/include/eli/util/tolerance.hpp
                                   ${ELI_SOURCE_DIR}/include/eli/geom/general/continuity.hpp
                                   ${ELI_SOURCE_DIR}/include/eli/geom/general/bounding_box.hpp
                                   ${ELI_SOURCE_DIR}/include/eli/geom/point/distance.hpp
                                   ${ELI_SOURCE_DIR}/include/eli/geom/curve/bezier.hpp
                                   ${ELI_SOURCE_DIR}/include/eli/geom/curve/equivalent_curves.hpp
                                   ${ELI_SOURCE_DIR}/include/eli/geom/surface/bezier.hpp
                                   ${ELI_SOURCE_DIR}/include/eli/geom/surface/curvature.hpp
                                   ${ELI_SOURCE_DIR}/include/eli/geom/surface/piecewise.hpp)

add_executable(${PIECEWISE_SURFACE_TEST_EXE} ${PIECEWISE_SURFACE_TEST_SOURCE} ${PIECEWISE_SURFACE_TEST_HEADER})
target_link_libraries(${PIECEWISE_SURFACE_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${PIECEWISE_SURFACE_TEST_EXE} COMMAND $<TARGET_FILE:${PIECEWISE_SURFACE_TEST_EXE}> --text)

# add the piecewise_surface_creator test
set (PIECEWISE_SURFACE_CREATOR_TEST_EXE PiecewiseSurfaceCreatorTest)
set (PIECEWISE_SURFACE_CREATOR_TEST_SOURCE ${PROJECT_SOURCE_DIR}/piecewise_surface_creator_test.cpp)
set (PIECEWISE_SURFACE_CREATOR_TEST_HEADER ${PROJECT_SOURCE_DIR}/piecewise_body_of_revolution_creator_test_suite.hpp
                                           ${PROJECT_SOURCE_DIR}/piecewise_skinning_surface_creator_test_suite.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/util/tolerance.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/general/continuity.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/general/bounding_box.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/point/distance.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/curve/bezier.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/curve/equivalent_curves.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/curve/piecewise.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/surface/bezier.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/surface/curvature.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/surface/piecewise.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/surface/piecewise_body_of_revolution_creator.hpp
                                           ${ELI_SOURCE_DIR}/include/eli/geom/surface/piecewise_skinning_surface_creator.hpp)

add_executable(${PIECEWISE_SURFACE_CREATOR_TEST_EXE} ${PIECEWISE_SURFACE_CREATOR_TEST_SOURCE} ${PIECEWISE_SURFACE_CREATOR_TEST_HEADER})
target_link_libraries(${PIECEWISE_SURFACE_CREATOR_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${PIECEWISE_SURFACE_CREATOR_TEST_EXE} COMMAND $<TARGET_FILE:${PIECEWISE_SURFACE_CREATOR_TEST_EXE}> --text)
