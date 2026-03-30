# (C) Copyright 2026 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.



# Set compiler flags for basic build types,
# for compilers where this is not provided by ecbuild.
include(build_type_compiler_flags)

# Set JEDI's common compiler flags
include(jedi_common_compiler_flags)

# Set UFO-specific compiler flags
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  ecbuild_add_fortran_flags("-ffree-line-length-none")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES Intel)  # Intel or IntelLLVM
  # some UFO code is sensitive to fp-model
  ecbuild_add_cxx_flags("-fp-model=strict")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)  # Intel or IntelLLVM
  # some UFO code is sensitive to fp-model
  ecbuild_add_fortran_flags("-fp-model=strict")
endif()
