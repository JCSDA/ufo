list( APPEND ufo_src_files
Fortran.h
GeoVaLs.cc
GeoVaLs.h
Locations.h
ObsBias.h
ObsSpace.cc
ObsSpace.h
ObsVector.cc
ObsVector.h
ObsWSpeed.cc
ObsWSpeed.h
ObsRadiance.cc
ObsRadiance.h
UfoTrait.h
Variables.h
Variables.cc
ufo_geovals_mod.F90
ufo_geovals_mod_c.F90
ufo_locs_mod.F90
ufo_obs_data.F90
ufo_obs_vectors.F90
ufo_variables_mod.F90
ufo_variables_mod_c.F90
ufo_wspeed_mod.F90
ufo_radiance_mod.F90
)

ecbuild_add_library( TARGET   ufo
                     SOURCES  ${ufo_src_files}
                     LIBS     oops ioda ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${CRTM_LIBRARIES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

