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

project(AIRFOIL_FIT_EXAMPLE)

# add the airfoil fitting example
set (AIRFOIL_FIT_EXE AirfoilFitExample)
set (AIRFOIL_FIT_SOURCE ${PROJECT_SOURCE_DIR}/airfoil_fit.cpp)
set (AIRFOIL_FIT_HEADER ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                        ${ELI_SOURCE_DIR}/include/eli/geom/curve/bezier.hpp)

add_executable(${AIRFOIL_FIT_EXE} ${AIRFOIL_FIT_SOURCE} ${AIRFOIL_FIT_HEADER})
target_link_libraries(${AIRFOIL_FIT_EXE} ${_ELI_EXAMPLE_LIBRARIES_LIST})
