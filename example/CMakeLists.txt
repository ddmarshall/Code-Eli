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

project(ELI_EXAMPLE)

# set what libraries need to be used
set(_ELI_EXAMPLE_LIBRARIES_LIST "")
if (QD_FOUND)
  SET (_ELI_EXAMPLE_LIBRARIES_LIST ${_ELI_EXAMPLE_LIBRARIES_LIST} ${QD_LIBRARIES})
endif()

add_subdirectory(airfoil_fit)
