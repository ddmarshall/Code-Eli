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
add_executable(VSPPodExample ${PROJECT_SOURCE_DIR}/vsp_pod.cpp)
target_link_libraries(VSPPodExample ${_ELI_EXAMPLE_LIBRARIES_LIST})
