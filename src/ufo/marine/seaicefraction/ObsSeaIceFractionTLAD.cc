/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/marine/seaicefraction/ObsSeaIceFractionTLAD.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsSeaIceFractionTLAD> makerSeaIceFractionTL_("SeaIceFraction");
// -----------------------------------------------------------------------------

ObsSeaIceFractionTLAD::ObsSeaIceFractionTLAD(const ioda::ObsSpace & odb,
                                             const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"sea_ice_category_area_fraction"};
  varin_.reset(new oops::Variables(vv));
  std::cout << keyOper_ << std::endl;
  ufo_seaicelinear_setup_f90(keyOper_, config);
  std::cout << keyOper_ << std::endl;
  oops::Log::trace() << "ObsSeaIceFractionTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFractionTLAD::~ObsSeaIceFractionTLAD() {
  ufo_seaicelinear_delete_f90(keyOper_);
  oops::Log::trace() << "ObsSeaIceFractionTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                          ObsDiagnostics &) {
  ufo_seaicelinear_settraj_f90(keyOper_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsSeaIceFractionTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::simulateObsTL(const GeoVaLs & gv, ioda::ObsVector & ovec) const {
  int nlocs = ovec.size();
  int nlevs = gv.nlevs("sea_ice_category_area_fraction");

  std::vector<double> aicen(nlocs);
  for ( std::size_t k = 1; k < nlevs+1; ++k ) {
    gv.get(aicen, "sea_ice_category_area_fraction", k);
    for ( std::size_t i = 0; i < nlocs; ++i ) {
      ovec[i] += aicen[i];
    }
  }
  oops::Log::trace() << "ObsSeaIceFractionTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::simulateObsAD(GeoVaLs & gv, const ioda::ObsVector & ovec) const {
  ufo_seaicelinear_alloc_ad_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran());
  int nlocs = ovec.size();
  int nlevs = gv.nlevs("sea_ice_category_area_fraction");
  float miss = 0.0;
  std::vector<double> aicen(nlocs);
  for ( std::size_t k = 1; k < nlevs+1; ++k ) {
    for ( std::size_t i = 0; i < nlocs; ++i ) {
      if (ovec[i] != util::missingValue(ovec[i]))
        { aicen[i] = ovec[i]; }
      else
        { aicen[i] = 0.0; }
    }
    gv.put(aicen, "sea_ice_category_area_fraction", k);
  }
  oops::Log::trace() << "ObsSeaIceFractionTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::print(std::ostream & os) const {
  os << "ObsSeaIceFractionTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
