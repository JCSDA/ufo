/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/marine/marinevertinterp/ObsMarineVertInterpTLAD.h"

#include <ostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsMarineVertInterpTLAD> makerMarinevertinterpTL_("InsituSalinity");
// -----------------------------------------------------------------------------

ObsMarineVertInterpTLAD::ObsMarineVertInterpTLAD(const ioda::ObsSpace & odb,
                                                   const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_(), varout_()
{
  int c_name_size = 200;
  char *buffin = new char[c_name_size];
  char *buffout = new char[c_name_size];
  const eckit::Configuration * configc = &config;
  ufo_marinevertinterp_tlad_setup_f90(keyOper_, &configc, buffin, buffout, c_name_size);

  std::string vstr_in(buffin), vstr_out(buffout);
  std::vector<std::string> vvin;
  std::vector<std::string> vvout;
  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
  varin_.reset(new oops::Variables(vvin));

  oops::Log::trace() << "ObsMarineVertInterpTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsMarineVertInterpTLAD::~ObsMarineVertInterpTLAD() {
  ufo_marinevertinterp_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsMarineVertInterpTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_marinevertinterp_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_);
  oops::Log::trace() << "ObsMarineVertInterpTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                             const ObsBiasIncrement & bias) const {
  ufo_marinevertinterp_simobs_tl_f90(keyOper_, geovals.toFortran(), odb_,
                                      ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsMarineVertInterpTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                             ObsBiasIncrement & bias) const {
  ufo_marinevertinterp_simobs_ad_f90(keyOper_, geovals.toFortran(), odb_,
                                      ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsMarineVertInterpTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsMarineVertInterpTLAD::print(std::ostream & os) const {
  os << "ObsMarinevertinterpTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
