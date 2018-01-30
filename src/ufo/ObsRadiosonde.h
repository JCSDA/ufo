/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSRADIOSONDE_H_
#define UFO_OBSRADIOSONDE_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsOperatorBase.h"
#include "ObsSpace.h"
#include "GeoVaLs.h"
#include "Locations.h"
#include "ObsBias.h"
#include "ObsBiasIncrement.h"
#include "ObsVector.h"
#include "util/ObjectCounter.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Radiosonde (currently only temperature) observation for UFO.
template <typename MODEL>
class ObsRadiosonde : public oops::ObsOperatorBase<MODEL>,
                  private util::ObjectCounter<ObsRadiosonde<MODEL>> {
 public:
  static const std::string classname() {return "ufo::ObsRadiosonde";}

  ObsRadiosonde(const ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadiosonde();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperRadiosonde_;}
  const int & toFortran() const {return keyOperRadiosonde_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperRadiosonde_;
  const ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsRadiosonde<MODEL>::ObsRadiosonde(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOperRadiosonde_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_radiosonde_setup_f90(keyOperRadiosonde_, &configc);
  const std::vector<std::string> vv{"virtual_temperature", "atmosphere_ln_pressure_coordinate"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsRadiosonde created." << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsRadiosonde<MODEL>::~ObsRadiosonde() {
  ufo_radiosonde_delete_f90(keyOperRadiosonde_);
  oops::Log::trace() << "ObsRadiosonde destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsRadiosonde<MODEL>::obsEquiv(const GeoVaLs & gom, ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_radiosonde_t_eqv_f90(gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsRadiosonde<MODEL>::print(std::ostream & os) const {
  os << "ObsRadiosonde::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSRADIOSONDE_H_
