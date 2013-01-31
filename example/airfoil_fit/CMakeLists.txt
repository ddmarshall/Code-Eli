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

project(AIRFOIL_FIT)

# add the airfoil fitting example
add_executable(AirfoilFitExample ${PROJECT_SOURCE_DIR}/airfoil_fit.cpp)
target_link_libraries(AirfoilFitExample ${_ELI_EXAMPLE_LIBRARIES_LIST})
