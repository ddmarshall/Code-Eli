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

project(VSP_POD)

# add the OpenVSP Pod construction example
set (VSP_POD_EXE VSPPodExample)
set (VSP_POD_SOURCE ${PROJECT_SOURCE_DIR}/vsp_pod.cpp)
set (VSP_POD_HEADER ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                    ${ELI_SOURCE_DIR}/include/eli/geom/curve/bezier.hpp
                    ${ELI_SOURCE_DIR}/include/eli/geom/curve/piecewise.hpp
                    ${ELI_SOURCE_DIR}/include/eli/geom/curve/piecewise_creator.hpp
                    ${ELI_SOURCE_DIR}/include/eli/geom/surface/bezier.hpp
                    ${ELI_SOURCE_DIR}/include/eli/geom/surface/piecewise.hpp
                    ${ELI_SOURCE_DIR}/include/eli/geom/surface/piecewise_creator.hpp)

add_executable(${VSP_POD_EXE} ${VSP_POD_SOURCE} ${VSP_POD_HEADER})
target_link_libraries(${VSP_POD_EXE} ${_ELI_EXAMPLE_LIBRARIES_LIST})
