set ( basis_files
	ufo_basis_mod.F90
	ufo_conventional_profile_mod.F90
	ufo_basis_tlad_mod.F90
	ufo_conventional_profile_tlad_mod.F90
)

PREPEND( _p_basis_files       "basis"       ${basis_files} )

set ( basis_src_files
    ${_p_basis_files} 
	PARENT_SCOPE
)
