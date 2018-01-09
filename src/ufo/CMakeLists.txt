list( APPEND ufo_src_files
Fortran.h
GeoVaLs.cc
GeoVaLs.h
Locations.cc
Locations.h
ObsBias.h
ObsBiasCovariance.h
ObsBiasIncrement.h
ObsSpace.cc
ObsSpace.h
ObsVector.cc
ObsVector.h
ObsSeaIceFraction.cc
ObsSeaIceFraction.h
UfoTrait.h
random_f.cc
random_f.h
random_vectors_mod.F90
GeoVaLs.interface.F90
ufo_geovals_mod.F90
Locations.interface.F90
ufo_locs_mod.F90
ufo_variables_mod.F90
ufo_obs_vectors.F90
ObsSeaIceFraction.interface.F90
ufo_seaicefrac_mod.F90
ObsSpace.SeaIceFraction.interface.F90
ufo_obs_seaicefrac_mod.F90
)

ecbuild_add_library( TARGET   ufo
                     SOURCES  ${ufo_src_files}
                     LIBS     oops ioda ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${CRTM_LIBRARIES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

