# Configures the compiler flags
#
# For each supported compiler vendor, the compiler settings for C, C++ and Fortran are
# set. This includes the setting of general flags, debug flags, release flags, minsize
# release flags and release with debug flags.
#
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

# only do this once
if(NOT CONFIGURE_COMPILER_INCLUDED)

  # set flag so will not process again
  set(CONFIGURE_COMPILER_INCLUDED TRUE)

  # flag to know if compiler flags were set
  set(CONFIGURE_COMPILER_SET_C_COMPILER_FLAGS OFF)
  set(CONFIGURE_COMPILER_SET_CXX_COMPILER_FLAGS OFF)
  set(CONFIGURE_COMPILER_SET_Fortran_COMPILER_FLAGS OFF)

  # set Intel compiler flags
  if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
    if(WIN32)
      set(CONFIGURE_COMPILER_SET_C_COMPILER_FLAGS ON)
      set(CMAKE_C_FLAGS "/Qstd=c99"
                        CACHE STRING "Flags used by the compiler during all build types." FORCE)
      set(CMAKE_C_FLAGS_DEBUG "/DDEBUG /Zi"
                        CACHE STRING "Flags used by the compiler during debug builds." FORCE)
      set(CMAKE_C_FLAGS_RELEASE "/DNDEBUG /O3 /ipo /fast"
                        CACHE STRING "Flags used by the compiler during release builds." FORCE)
      set(CMAKE_C_FLAGS_RELWITHDEBINFO "/DDEBUG /O3 /Zi"
                        CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
      set(CMAKE_C_FLAGS_MINSIZEREL "/DNDEBUG /Os /ipo /fast"
                        CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
    else()
      set(CONFIGURE_COMPILER_SET_C_COMPILER_FLAGS ON)
      set(CMAKE_C_FLAGS "-std=c99"
                        CACHE STRING "Flags used by the compiler during all build types." FORCE)
      set(CMAKE_C_FLAGS_DEBUG "-DDEBUG -O0 -g"
                        CACHE STRING "Flags used by the compiler during debug builds." FORCE)
      set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -O3 -ipo -fast -g0"
                        CACHE STRING "Flags used by the compiler during release builds." FORCE)
      set(CMAKE_C_FLAGS_RELWITHDEBINFO "-DDEBUG -O3 -ipo -g"
                        CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
      set(CMAKE_C_FLAGS_MINSIZEREL "-DNDEBUG -Os -ipo -fast"
                        CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)

      # Hack because CMake does not set this flag for intel fortran compiler
      set(CMAKE_Fortran_COMPILER_ID "Intel")
    endif()
  endif()
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    if(WIN32)
      set(CONFIGURE_COMPILER_SET_CXX_COMPILER_FLAGS ON)
      set(CMAKE_CXX_FLAGS "/Qstd=c++0x"
                          CACHE STRING "Flags used by the compiler during all build types." FORCE)
      set(CMAKE_CXX_FLAGS_DEBUG "/DDEBUG /Zi"
                          CACHE STRING "Flags used by the compiler during debug builds." FORCE)
      set(CMAKE_CXX_FLAGS_RELEASE "/DNDEBUG /DEIGEN_NO_DEBUG /O3 /ipo /fast"
                          CACHE STRING "Flags used by the compiler during release builds." FORCE)
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/DDEBUG /O3 /Zi"
                          CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
      set(CMAKE_CXX_FLAGS_MINSIZEREL "/DNDEBUG /DEIGEN_NO_DEBUG /Os /ipo /fast"
                          CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
    else()
      set(CONFIGURE_COMPILER_SET_CXX_COMPILER_FLAGS ON)
      set(CMAKE_CXX_FLAGS "-std=c++0x"
                          CACHE STRING "Flags used by the compiler during all build types." FORCE)
      set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -O0 -g"
                          CACHE STRING "Flags used by the compiler during debug builds." FORCE)
      set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -DEIGEN_NO_DEBUG -O3 -ipo -fast -g0"
                          CACHE STRING "Flags used by the compiler during release builds." FORCE)
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DDEBUG -O3 -ipo -g"
                          CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
      set(CMAKE_CXX_FLAGS_MINSIZEREL "-DNDEBUG -DEIGEN_NO_DEBUG -Os -ipo -fast"
                          CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
    endif()
  endif()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if(WIN32)
      set(CONFIGURE_COMPILER_SET_Fortran_COMPILER_FLAGS OFF)
    else()
      set(CONFIGURE_COMPILER_SET_Fortran_COMPILER_FLAGS ON)
      set(CMAKE_Fortran_FLAGS ""
                              CACHE STRING "Flags used by the compiler during all build types." FORCE)
      set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g"
                              CACHE STRING "Flags used by the compiler during debug builds." FORCE)
      set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo"
                              CACHE STRING "Flags used by the compiler during release builds." FORCE)
      set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -ipo -g"
                              CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
      set(CMAKE_Fortran_FLAGS_MINSIZEREL "-Os -ipo"
                              CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
    endif()
  endif()

  # set MS Visual 2008 Compiler flags
  if (MSVC90)
    if (CMAKE_CL_64)
      set(CONFIGURE_COMPILER_SET_C_COMPILER_FLAGS OFF)
      set(CONFIGURE_COMPILER_SET_CXX_COMPILER_FLAGS OFF)
    else()
      set(CONFIGURE_COMPILER_SET_C_COMPILER_FLAGS ON)
      set(CMAKE_C_FLAGS "/DWIN32 /D_WINDOWS /W3 /Zm1000 /Qc99"
                        CACHE STRING "Flags used by the compiler during all build types." FORCE)
      set(CMAKE_C_FLAGS_DEBUG "/D_DEBUG /MDd /Zi  /Ob0 /Od /RTC1 /DDEBUG"
                        CACHE STRING "Flags used by the compiler during debug builds." FORCE)
      set(CMAKE_C_FLAGS_RELEASE "/MD /O2 /Ob2 /DNDEBUG"
                        CACHE STRING "Flags used by the compiler during release builds." FORCE)
      set(CMAKE_C_FLAGS_RELWITHDEBINFO "/MD /Zi /O2 /Ob1 /DNDEBUG"
                        CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
      set(CMAKE_C_FLAGS_MINSIZEREL "/MD /O1 /Ob1 /DNDEBUG"
                        CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)

      set(CONFIGURE_COMPILER_SET_CXX_COMPILER_FLAGS ON)
      set(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS /W3 /Zm1000 /EHsc /GR"
                          CACHE STRING "Flags used by the compiler during all build types." FORCE)
      set(CMAKE_CXX_FLAGS_DEBUG "/D_DEBUG /MDd /Zi /Ob0 /Od /RTC1 /DDEBUG"
                          CACHE STRING "Flags used by the compiler during debug builds." FORCE)
      set(CMAKE_CXX_FLAGS_RELEASE "/MD /O2 /Ob2 /DNDEBUG /DEIGEN_NO_DEBUG"
                          CACHE STRING "Flags used by the compiler during release builds." FORCE)
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MD /Zi /O2 /Ob1 /DNDEBUG"
                          CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
      set(CMAKE_CXX_FLAGS_MINSIZEREL "/MD /O1 /Ob1 /D NDEBUG /DEIGEN_NO_DEBUG"
                          CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
    endif()

    set(CONFIGURE_COMPILER_SET_Fortran_COMPILER_FLAGS OFF)
  endif()

  # set GNU Compiler Collection flags
  if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set(CONFIGURE_COMPILER_SET_C_COMPILER_FLAGS ON)
    set(CMAKE_C_FLAGS "-ansi -pedantic -Wall -Wextra -fmessage-length=100 -std=c99"
                      CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_C_FLAGS_DEBUG "-DDEBUG -O0 -g"
                      CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -O3"
                      CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-DDEBUG -O3 -g"
                      CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
    set(CMAKE_C_FLAGS_MINSIZEREL "-DNDEBUG -Os"
                      CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)

    # Hack because CMake does not set this flag for GNU fortran compiler
    set(CMAKE_Fortran_COMPILER_ID "GNU")
  endif()
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CONFIGURE_COMPILER_SET_CXX_COMPILER_FLAGS ON)
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
      if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
        set(CMAKE_CXX_FLAGS "-std=c++11")
      else()
        set(CMAKE_CXX_FLAGS "-std=c++0x")
      endif()
    else()
      message(FATAL_ERROR "Version ${GCC_VERSION} of GCC compiler is not supported.")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -pedantic -Wall -Wextra -Wno-long-long -fmessage-length=100"
                        CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -O0 -g"
                        CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -DEIGEN_NO_DEBUG -O3"
                        CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DDEBUG -O3 -g"
                        CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
    set(CMAKE_CXX_FLAGS_MINSIZEREL "-DNDEBUG -DEIGEN_NO_DEBUG -Os"
                        CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
  endif()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CONFIGURE_COMPILER_SET_Fortran_COMPILER_FLAGS ON)
    set(CMAKE_Fortran_FLAGS ""
                            CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_Fortran_FLAGS_DEBUG "-DDEBUG -O0 -g"
                            CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG -O3"
                            CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-DDEBUG -O3 -g"
                            CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
    set(CMAKE_Fortran_FLAGS_MINSIZEREL "-DNDEBUG -Os"
                            CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
  endif()

  # set LLVM flags
  if(CMAKE_C_COMPILER_ID STREQUAL "Clang")
    set(CONFIGURE_COMPILER_SET_C_COMPILER_FLAGS ON)
    if (CMAKE_GENERATOR STREQUAL "Xcode")
      set(CMAKE_XCODE_ATTRIBUTE_GCC_C_LANGUAGE_STANDARD "c99")
      set(CMAKE_XCODE_ATTRIBUTE_GCC_WARN_INHIBIT_ALL_WARNINGS "YES")
      set(CMAKE_XCODE_ATTRIBUTE_GCC_WARN_PEDANTIC "YES")
    else()
      set(CMAKE_C_FLAGS "-pedantic -Wall -Wextra -fmessage-length=100 -std=c99"
                        CACHE STRING "Flags used by the compiler during all build types." FORCE)
    endif()
    set(CMAKE_C_FLAGS_DEBUG "-DDEBUG -O0 -g"
                      CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -O3"
                      CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-DDEBUG -O3 -g"
                      CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
    set(CMAKE_C_FLAGS_MINSIZEREL "-DNDEBUG -Os"
                      CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)

    # Hack because CMake does not set this flag for Clang fortran compiler
    set(CMAKE_Fortran_COMPILER_ID "Clang")
  endif()
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CONFIGURE_COMPILER_SET_CXX_COMPILER_FLAGS ON)
    if (CMAKE_GENERATOR STREQUAL "Xcode")
      set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
      set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
      set(CMAKE_XCODE_ATTRIBUTE_GCC_WARN_INHIBIT_ALL_WARNINGS "YES")
      set(CMAKE_XCODE_ATTRIBUTE_GCC_WARN_PEDANTIC "YES")
    else()
      set(CMAKE_CXX_FLAGS "-pedantic -Wall -Wextra -Wno-long-long -fmessage-length=100 -std=c++11 -stdlib=libc++"
                          CACHE STRING "Flags used by the compiler during all build types." FORCE)
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -O0 -g"
                        CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    # FIX: SET SOURCE FILE PROPERTIES ON AFFECTED FILES INSTEAD OF THIS
    #      (1) Try and figure out problem and fix code
    #      (2) SET_SOURCE_FILES_PROPERTIES(${file} PROPERTIES COMPILE_FLAGS -O1)
    #      Clang 3.3 and lower produces wrong answer for mutilsimpsontestsuite hangs for -O2 or -O3
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -DEIGEN_NO_DEBUG -O3"
                        CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DDEBUG -O3 -g"
                        CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
      set(CMAKE_CXX_FLAGS_MINSIZEREL "-DNDEBUG -DEIGEN_NO_DEBUG -Os"
                        CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
  endif()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Clang")
    set(CONFIGURE_COMPILER_SET_Fortran_COMPILER_FLAGS ON)
    set(CMAKE_Fortran_FLAGS                ""
                      CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_Fortran_FLAGS_DEBUG "-DDEBUG -O0 -g"
                      CACHE STRING "Flags used by the compiler during debug builds." FORCE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG -O3"
                      CACHE STRING "Flags used by the compiler during release builds." FORCE)
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-DDEBUG -O3 -g"
                      CACHE STRING "Flags used by the compiler during release with debug info builds." FORCE)
      set(CMAKE_Fortran_FLAGS_MINSIZEREL "/DNDEBUG /Os /ipo /fast"
                      CACHE STRING "Flags used by the compiler during release minsize builds." FORCE)
  endif()

  if (NOT CONFIGURE_COMPILER_SET_C_COMPILER_FLAGS)
    message(FATAL_ERROR "Could not set C Compiler flags! \nNeed to fix ConfigureCompiler.cmake to support ${CMAKE_C_COMPILER_ID} compiler!")
  endif()
  if (NOT CONFIGURE_COMPILER_SET_CXX_COMPILER_FLAGS)
    message(FATAL_ERROR "Could not set C++ Compiler flags! \nNeed to fix ConfigureCompiler.cmake to support ${CMAKE_CXX_COMPILER_ID} compiler!")
  endif()
#  if (NOT CONFIGURE_COMPILER_SET_Fortran_COMPILER_FLAGS)
#    message(FATAL_ERROR "Could not set Fortran Compiler flags! \nNeed to fix ConfigureCompiler.cmake to support ${CMAKE_Fortran_COMPILER_ID} compiler!")
#  endif()
endif()
