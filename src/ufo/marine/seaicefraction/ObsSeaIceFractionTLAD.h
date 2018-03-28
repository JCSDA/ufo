/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSSEAICEFRACTIONTLAD_H_
#define UFO_OBSSEAICEFRACTIONTLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperBase.h"
#include "ufo/ObsSpace.h"
#include "util/ObjectCounter.h"
#include "util/Logger.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVector;

// -----------------------------------------------------------------------------
/// Sea-ice fraction observation for  model.
template <typename MODEL>
class ObsSeaIceFractionTLAD : public oops::LinearObsOperBase<MODEL>, 
                              private util::ObjectCounter<ObsSeaIceFractionTLAD<MODEL>> {
public:
  static const std::string classname() {return "ufo::ObsSeaIceFractionTLAD";}

  ObsSeaIceFractionTLAD(const ObsSpace &, const eckit::Configuration &);    
  virtual ~ObsSeaIceFractionTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void obsEquivTL(const GeoVaLs &, ObsVector &, const ObsBiasIncrement &) const;
  void obsEquivAD(GeoVaLs &, const ObsVector &, ObsBiasIncrement &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperSeaIceFraction_;}
  const int & toFortran() const {return keyOperSeaIceFraction_;}

private:
  void print(std::ostream &) const;
  F90hop keyOperSeaIceFraction_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsSeaIceFractionTLAD<MODEL>::ObsSeaIceFractionTLAD(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperSeaIceFraction_(0), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_seaicefrac_setup_f90(keyOperSeaIceFraction_, &configc);
  const std::vector<std::string> vv{"ice_concentration"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsSeaIceFractionTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsSeaIceFractionTLAD<MODEL>::~ObsSeaIceFractionTLAD() {
  oops::Log::trace() << "ObsSeaIceFractionTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsSeaIceFractionTLAD<MODEL>::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_seaicefrac_settraj_f90(keyOperSeaIceFraction_, geovals.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsSeaIceFractionTLAD<MODEL>::obsEquivTL(const GeoVaLs & geovals, ObsVector & ovec,
                               const ObsBiasIncrement & bias) const {
  ufo_seaicefrac_eqv_tl_f90(keyOperSeaIceFraction_, geovals.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsSeaIceFractionTLAD<MODEL>::obsEquivAD(GeoVaLs & geovals, const ObsVector & ovec,
                               ObsBiasIncrement & bias) const {
  ufo_seaicefrac_eqv_ad_f90(keyOperSeaIceFraction_, geovals.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsSeaIceFractionTLAD<MODEL>::print(std::ostream & os) const {
  os << "ObsSeaIceFractionTLAD::print not implemented" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSSEAICEFRACTIONTLAD_H_
