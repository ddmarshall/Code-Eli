# Determine the date and time when called. Can be called multiple times to get the
# date and time at multiple points in CMake process.
#
# This tries to find the CPPTest library. Once completed it will set:
# ELI_DATE - date in YYYYMMDD format
# ELI_TIME - time in HHMMSS format
# ELI_DATE_TIME_FOUND - boolean set to true if was able to determine date and time
#
################################################################################
# Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
#
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the Eclipse Public License v1.0
# which accompanies this distribution, and is available at
# http://www.eclipse.org/legal/epl-v10.html
#
# Contributors:
#    David D. Marshall - initial code and implementation
################################################################################

# reset values in case already run.
set(ELI_DATE "DATE_NOT_FOUND")
set(ELI_TIME "TIME_NOT_FOUND")
set(ELI_DATE_TIME_FOUND FALSE)

if (WIN32)
  # Calculate the date
  execute_process(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE _ELI_WIN_DATE)
  string(REGEX MATCH "(..)/(..)/(....)" _ELI_WIN_DATE ${_ELI_WIN_DATE})
  string(REGEX REPLACE "(..)/(..)/(....)" "\\3\\1\\2" ELI_DATE ${_ELI_WIN_DATE})

  # Calculate the time
  execute_process(COMMAND "cmd" " /C echo %time%" OUTPUT_VARIABLE _ELI_WIN_TIME)
  string(REGEX MATCH "(..):(..):(..)" _ELI_WIN_TIME ${_ELI_WIN_TIME})
  string(REGEX REPLACE "(..):(..):(..)" "\\1" _ELI_WIN_HOUR_TIME ${_ELI_WIN_TIME})
  string(STRIP ${_ELI_WIN_HOUR_TIME} _ELI_WIN_HOUR_TIME)
  if(_ELI_WIN_HOUR_TIME LESS "10")
    set(_ELI_WIN_HOUR_TIME "0${_ELI_WIN_HOUR_TIME}")
  endif()
  string(REGEX REPLACE "(..):(..):(..)" "\\2\\3" _ELI_WIN_MINSEC_TIME ${_ELI_WIN_TIME})
  set(ELI_TIME "${_ELI_WIN_HOUR_TIME}${_ELI_WIN_MINSEC_TIME}")
  set(ELI_DATE_TIME_FOUND TRUE)
elseif (UNIX)
  # Calculate the date
  execute_process(COMMAND "date" "+%Y%m%d" OUTPUT_VARIABLE ELI_UNIX_DATE)
  string(REGEX REPLACE "(\r?\n)+$" "" ELI_DATE "${ELI_UNIX_DATE}")

  # Calculate the time
  execute_process(COMMAND "date" "+%H%M%S" OUTPUT_VARIABLE ELI_UNIX_TIME)
  string(REGEX REPLACE "(\r?\n)+$" "" ELI_TIME "${ELI_UNIX_TIME}")
  set(ELI_DATE_TIME_FOUND TRUE)
else()
endif()
