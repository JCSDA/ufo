/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSRADIOSONDETLAD_H_
#define UFO_OBSRADIOSONDETLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperBase.h"
#include "ioda/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Logger.h"
#include "ufo/FortranAtmosphere.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// Radiosonde (currently only temperature) observation for UFO.
template <typename MODEL>
class ObsRadiosondeTLAD : public oops::LinearObsOperBase<MODEL>,
                          private util::ObjectCounter<ObsRadiosondeTLAD<MODEL>> {
public:
  static const std::string classname() {return "ufo::ObsRadiosondeTLAD";}

  ObsRadiosondeTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadiosondeTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void obsEquivTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void obsEquivAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperRadiosonde_;}
  const int & toFortran() const {return keyOperRadiosonde_;}

private:
  void print(std::ostream &) const;
  F90hop keyOperRadiosonde_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsRadiosondeTLAD<MODEL>::ObsRadiosondeTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperRadiosonde_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_radiosonde_tlad_setup_f90(keyOperRadiosonde_, &configc);
  const std::vector<std::string> vv{"virtual_temperature", "atmosphere_ln_pressure_coordinate"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsRadiosondeTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsRadiosondeTLAD<MODEL>::~ObsRadiosondeTLAD() {
  ufo_radiosonde_tlad_delete_f90(keyOperRadiosonde_);
  oops::Log::trace() << "ObsRadiosondeTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsRadiosondeTLAD<MODEL>::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_radiosonde_tlad_settraj_f90(keyOperRadiosonde_, geovals.toFortran(), odb_.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsRadiosondeTLAD<MODEL>::obsEquivTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                             const ObsBiasIncrement & bias) const {
  ufo_radiosonde_tlad_t_eqv_tl_f90(keyOperRadiosonde_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsRadiosondeTLAD<MODEL>::obsEquivAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                             ObsBiasIncrement & bias) const {
  ufo_radiosonde_tlad_t_eqv_ad_f90(keyOperRadiosonde_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsRadiosondeTLAD<MODEL>::print(std::ostream & os) const {
  os << "ObsRadiosondeTLAD::print not implemented" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSRADIOSONDETLAD_H_
