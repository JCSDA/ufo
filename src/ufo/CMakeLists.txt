list( APPEND ufo_src_files
    Fortran.h
    GeoVaLs.cc
    GeoVaLs.h
    ObsBias.h
    ObsBiasCovariance.h
    ObsBiasIncrement.h
    ObsCheck.cc
    ObsCheck.h
    UfoTrait.h
    instantiateObsOperatorFactory.h
    instantiateLinearObsOpFactory.h
    ufo_obscheck_mod.F90
    GeoVaLs.interface.F90
    ufo_geovals_mod.F90
    ufo_variables_mod.F90
    ufo_obs_data_basis_mod.F90
    ufo_obs_data.F90
    ufo_obs_data_mod.F90
)

# macro to prepend a prefix with relative path
# can this be added to ecbuild for use elsewhere?
function(PREPEND var prefix )
    set ( listVar "" )
    foreach (f ${ARGN})
        list (APPEND listVar "${prefix}/${f}")
    endforeach(f)
    set ( ${var} "${listVar}" PARENT_SCOPE )
endfunction(PREPEND)

add_subdirectory( atmosphere )
add_subdirectory( constituents )
add_subdirectory( marine )

list( APPEND ufo_src_files
    ${atmosphere_src_files}
    ${constituents_src_files}
    ${marine_src_files}
)

message(STATUS ${ufo_src_files})
ecbuild_add_library( TARGET   ufo
                     SOURCES  ${ufo_src_files}
                     LIBS     oops ioda crtm gsw ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )
