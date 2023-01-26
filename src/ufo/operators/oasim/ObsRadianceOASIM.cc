/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/oasim/ObsRadianceOASIM.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"


namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceOASIM> makerOASIM_("OASIM");
// -----------------------------------------------------------------------------

ObsRadianceOASIM::ObsRadianceOASIM(const ioda::ObsSpace & odb, const Parameters_ & params)
  : ObsOperatorBase(odb), keyOper_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vvin{"surface_pressure_at_mean_sea_level",
                                      "surface_wind_speed",
                                      "ozone_thickness",
                                      "water_vapor",
                                      "relative_humidity",
                                      "cloud_area_fraction_in_atmosphere_layer",
                                      "cloud_optical_thickness",
                                      "cloud_liquid_water_path",
                                      "effective_radius_of_cloud_liquid_water_particle",
                                      "aerosol_optical_thickness",
                                      "single_scattering_albedo",
                                      "asymmetry_parameter",
                                      "sea_water_cell_thickness",
                                      "Carbon_nitrogen_detritus_concentration",
                                      "Particulate_inorganic_carbon",
                                      "colored_dissolved_organic_carbon",
                                      "diatom_concentration",
                                      "chlorophyte_concentration",
                                      "cyano-bacteria_concentration",
                                      "coccolithophore_concentration",
                                      "dinoflagellate_concentration",
                                      "phaeocystis_concentration",
                                      "cyano-bacteria_concentration"};

  varin_.reset(new oops::Variables(vvin));

  ufo_radianceoasim_setup_f90(keyOper_, params.toConfiguration());
  oops::Log::trace() << "ObsRadianceOASIM created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceOASIM::~ObsRadianceOASIM() {
  ufo_radianceoasim_delete_f90(keyOper_);
  oops::Log::trace() << "ObsRadianceOASIM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceOASIM::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              ObsDiagnostics &) const {
  ufo_radianceoasim_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsRadianceOASIM: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceOASIM::print(std::ostream & os) const {
  os << "ObsRadianceOASIM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
