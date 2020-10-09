# - Try to find SEP
# Once done this will define
#
#  SEP_FOUND - system has SEP
#  SEP_INCLUDE_DIR - the SEP include directory
#  SEP_LIBRARIES - Link these to use SEP
#  SEP_VERSION_STRING - Human readable version number of sep
#  SEP_VERSION_MAJOR  - Major version number of sep
#  SEP_VERSION_MINOR  - Minor version number of sep

# Copyright (c) 2017, Ilia Platone, <info@iliaplatone.com>
# Based on FindLibfacile by Carsten Niehaus, <cniehaus@gmx.de>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

if (SEP_LIBRARIES)

  # in cache already
  set(SEP_FOUND TRUE)
  message(STATUS "Found SEP: ${SEP_LIBRARIES}")


else (SEP_LIBRARIES)

  find_library(SEP_LIBRARIES NAMES sep
    PATHS
    ${_obLinkDir}
    ${GNUWIN32_DIR}/lib
    /usr/local/lib
  )

  if(SEP_LIBRARIES)
    set(SEP_FOUND TRUE)
  else (SEP_LIBRARIES)
    set(SEP_FOUND FALSE)
  endif(SEP_LIBRARIES)


  if (SEP_FOUND)
    if (NOT SEP_FIND_QUIETLY)
      message(STATUS "Found SEP: ${SEP_LIBRARIES}")
    endif (NOT SEP_FIND_QUIETLY)
  else (SEP_FOUND)
    if (SEP_FIND_REQUIRED)
      message(FATAL_ERROR "SEP not found. Please install libsep-dev")
    endif (SEP_FIND_REQUIRED)
  endif (SEP_FOUND)

  mark_as_advanced(SEP_LIBRARIES)
  
endif (SEP_LIBRARIES)
