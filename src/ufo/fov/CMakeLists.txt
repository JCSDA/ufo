# (C) Copyright 2024 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( fov_files
  FieldOfView.interface.F90
  FieldOfView.interface.h
  ReduceOverFieldOfView.h
  ReduceOverFieldOfView.cc
  SampleFieldOfView.h
  SampleFieldOfView.cc
  calc_fov_conical.f90
  calc_fov_crosstrk.f90
  ufo_fov_mod.f90
)

PREPEND( _p_fov_files "fov" ${fov_files} )

set ( fov_src_files
      ${_p_fov_files}
      PARENT_SCOPE
)

