# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( basis_files
	ufo_basis_mod.F90
	ufo_basis_tlad_mod.F90
)

PREPEND( _p_basis_files       "basis"       ${basis_files} )

set ( basis_src_files
    ${_p_basis_files} 
	PARENT_SCOPE
)
