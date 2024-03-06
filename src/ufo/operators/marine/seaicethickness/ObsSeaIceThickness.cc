/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/marine/seaicethickness/ObsSeaIceThickness.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsSeaIceThickness> makerSeaIceThickness_("SeaIceThickness");
// -----------------------------------------------------------------------------

ObsSeaIceThickness::ObsSeaIceThickness(const ioda::ObsSpace & odb,
                                       const ObsSeaIceThicknessParameters & params)
  : ObsOperatorBase(odb), keyOper_(0), odb_(odb), varin_()
{
  std::vector<std::string> vvin{"sea_ice_category_area_fraction",
                                "sea_ice_category_thickness"};
  if (odb.assimvariables().has("seaIceFreeboard")) {
    vvin.push_back("sea_ice_category_snow_thickness");
  }
  varin_.reset(new oops::Variables(vvin));
  ufo_seaicethickness_setup_f90(keyOper_, params.toConfiguration(), odb.assimvariables());
  oops::Log::trace() << "ObsSeaIceThickness created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceThickness::~ObsSeaIceThickness() {
  ufo_seaicethickness_delete_f90(keyOper_);
  oops::Log::trace() << "ObsSeaIceThickness destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceThickness::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                     ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  ufo_seaicethickness_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsSeaIceThickness: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceThickness::print(std::ostream & os) const {
  os << "ObsSeaIceThickness::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
