# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( ncdiag_files

ncd_kinds.F90
m_distribution.f90
m_diag_conv.f90
m_diag_aircft.f90
m_diag_raob.f90
nc_diag_fson.f90
nc_diag_read_mod.F90
nc_diag_res.f90
nc_diag_write_mod.F90
ncdc_cli_process.F90
ncdc_climsg.F90
ncdc_data.F90
ncdc_dims.F90
ncdc_metadata.F90
ncdc_realloc.F90
ncdc_state.F90
ncdc_types.f90
ncdc_vars.F90
ncdf_path_m.f90
ncdf_string_m.f90
ncdf_value_m.f90
ncdr_alloc_assert.f90
ncdr_attrs.f90
ncdr_attrs_fetch.f90
ncdr_check.f90
ncdr_climsg.F90
ncdr_dims.f90
ncdr_global_attrs.f90
ncdr_global_attrs_fetch.f90
ncdr_realloc_mod.F90
ncdr_state.f90
ncdr_types.f90
ncdr_vars.f90
ncdr_vars_fetch.f90
ncdres_climsg.F90
ncdw_chaninfo.F90
ncdw_ciresize.F90
ncdw_climsg.F90
ncdw_data2d.F90
ncdw_dresize.F90
ncdw_lheader.F90
ncdw_metadata.F90
ncdw_mresize.F90
ncdw_realloc.F90
ncdw_state.f90
ncdw_strarrutils.F90
ncdw_types.F90
ncdw_varattr.F90
netcdf_unlimdims.F90
read_diag.f90

)

PREPEND( _p_ncdiag_files       "ncdiag"       ${ncdiag_files} )

set ( ncdiag_src_files
      ${_p_ncdiag_files}
      PARENT_SCOPE
)
