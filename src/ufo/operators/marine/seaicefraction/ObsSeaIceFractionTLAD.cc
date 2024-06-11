/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/seaicefraction/ObsSeaIceFractionTLAD.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsSeaIceFractionTLAD> makerSeaIceFractionTL_("SeaIceFraction");
// -----------------------------------------------------------------------------

ObsSeaIceFractionTLAD::ObsSeaIceFractionTLAD(const ioda::ObsSpace & odb,
                                             const Parameters_ & params)
  : LinearObsOperatorBase(odb), varin_()
{
  const std::vector<std::string> vv{"sea_ice_category_area_fraction"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsSeaIceFractionTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFractionTLAD::~ObsSeaIceFractionTLAD() {
  oops::Log::trace() << "ObsSeaIceFractionTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                          const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsSeaIceFractionTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::simulateObsTL(const GeoVaLs & gv, ioda::ObsVector & ovec,
  const QCFlags_t& qc_flags) const {
  size_t nlocs = ovec.size();
  size_t nlevs = gv.nlevs(oops::Variable{"sea_ice_category_area_fraction"});

  std::vector<double> aicen(nlocs);
  for ( std::size_t k = 0; k < nlevs; ++k ) {
    gv.getAtLevel(aicen, oops::Variable{"sea_ice_category_area_fraction"}, k);
    for ( std::size_t i = 0; i < nlocs; ++i ) {
      ovec[i] += aicen[i];
    }
  }
  oops::Log::trace() << "ObsSeaIceFractionTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::simulateObsAD(GeoVaLs & gv, const ioda::ObsVector & ovec,
                                          const QCFlags_t & qc_flags) const {
  size_t nlocs = ovec.size();
  size_t nlevs = gv.nlevs(oops::Variable{"sea_ice_category_area_fraction"});
  std::vector<double> aicen(nlocs);
  for ( std::size_t k = 0; k < nlevs; ++k ) {
    for ( std::size_t i = 0; i < nlocs; ++i ) {
      if (ovec[i] != util::missingValue<double>())
        { aicen[i] = ovec[i]; }
      else
        { aicen[i] = 0.0; }
    }
    gv.putAtLevel(aicen, oops::Variable{"sea_ice_category_area_fraction"}, k);
  }
  oops::Log::trace() << "ObsSeaIceFractionTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFractionTLAD::print(std::ostream & os) const {
  os << "ObsSeaIceFractionTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
