/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/crtm/ObsAodCRTM.h"

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAodCRTM> makerAOD_("Aod");

// -----------------------------------------------------------------------------

ObsAodCRTM::ObsAodCRTM(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperAodCRTM_(0), odb_(odb), varin_(),
    channels_(odb.obsvariables().channels())
{
  const std::vector<std::string> vv{
    "air_temperature", "humidity_mixing_ratio", "relative_humidity",
    "air_pressure", "air_pressure_levels",
     "sulf", "bc1", "bc2", "oc1", "oc2", "dust1", "dust2", "dust3", "dust4", "dust5",
      "seas1", "seas2", "seas3", "seas4"};
  varin_.reset(new oops::Variables(vv));

  // call Fortran setup routine
  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");
  const eckit::Configuration * configc = &obsOptions;
  ufo_aodcrtm_setup_f90(keyOperAodCRTM_, &configc);
  oops::Log::info() << "ObsAodCRTM channels: " << channels_ << std::endl;
  oops::Log::trace() << "ObsAodCRTM created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodCRTM::~ObsAodCRTM() {
  ufo_aodcrtm_delete_f90(keyOperAodCRTM_);
  oops::Log::trace() << "ObsAodCRTM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCRTM::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec) const {
  ufo_aodcrtm_simobs_f90(keyOperAodCRTM_, gom.toFortran(), odb_,
                          ovec.size(), ovec.toFortran(),
                          channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsAodCRTM::print(std::ostream & os) const {
  os << "ObsAodCRTM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
