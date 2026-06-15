# (C) Copyright 2024-2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Author: Rong Kong (NCAR/MMM) - C++ rewrite of the PPRO operator from Fortran.

set ( ppro_src_files
    ObsPPRO.h
    ObsPPRO.cc
    PPROCoefsData.h
    Zhang21Forward.h
    Zhang21Forward.cc
    TCWA2Forward.h
    TCWA2Forward.cc
PARENT_SCOPE
)
