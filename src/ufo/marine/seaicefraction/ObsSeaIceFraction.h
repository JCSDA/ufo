/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSSEAICEFRACTION_H_
#define UFO_OBSSEAICEFRACTION_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsOperatorBase.h"
#include "ioda/ObsSpace.h"
#include "ufo/GeoVaLs.h"
#include "ioda/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ioda/ObsVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/FortranMarine.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Total ice concentration observation for UFO.
template <typename MODEL>
class ObsSeaIceFraction : public oops::ObsOperatorBase<MODEL>,
                  private util::ObjectCounter<ObsSeaIceFraction<MODEL>> {
 public:
  static const std::string classname() {return "ufo::ObsSeaIceFraction";}

  ObsSeaIceFraction(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsSeaIceFraction();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperSeaIceFraction_;}
  const int & toFortran() const {return keyOperSeaIceFraction_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperSeaIceFraction_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsSeaIceFraction<MODEL>::ObsSeaIceFraction(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperSeaIceFraction_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_seaicefrac_setup_f90(keyOperSeaIceFraction_, &configc);
  const std::vector<std::string> vv{"ice_concentration"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsSeaIceFraction created." << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsSeaIceFraction<MODEL>::~ObsSeaIceFraction() {
  ufo_seaicefrac_delete_f90(keyOperSeaIceFraction_);
  oops::Log::trace() << "ObsSeaIceFraction destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsSeaIceFraction<MODEL>::obsEquiv(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_seaicefrac_eqv_f90(keyOperSeaIceFraction_, gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsSeaIceFraction<MODEL>::print(std::ostream & os) const {
  os << "ObsSeaIceFraction::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSSEAICEFRACTION_H_
