/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/seaicefraction/ObsSeaIceFraction.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsSeaIceFraction> makerSeaIceFraction_("SeaIceFraction");
// -----------------------------------------------------------------------------

ObsSeaIceFraction::ObsSeaIceFraction(const ioda::ObsSpace & odb,
                                     const ObsSeaIceFractionParameters & params)
  : ObsOperatorBase(odb), varin_()
{
  const std::vector<std::string> vvin{"sea_ice_category_area_fraction"};
  varin_.reset(new oops::Variables(vvin));
  oops::Log::trace() << "ObsSeaIceFraction created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFraction::~ObsSeaIceFraction() {
  oops::Log::trace() << "ObsSeaIceFraction destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFraction::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                    ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  size_t nlocs = ovec.size();
  size_t nlevs = gv.nlevs(oops::Variable{"sea_ice_category_area_fraction"});

  std::vector<double> aicen(nlocs);
  for ( std::size_t k = 0; k < nlevs; ++k ) {
    gv.getAtLevel(aicen, oops::Variable{"sea_ice_category_area_fraction"}, k);
    for ( std::size_t i = 0; i < nlocs; ++i ) {
      ovec[i] += aicen[i];
    }
  }
  oops::Log::trace() << "ObsSeaIceFraction: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFraction::print(std::ostream & os) const {
  os << "ObsSeaIceFraction::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
