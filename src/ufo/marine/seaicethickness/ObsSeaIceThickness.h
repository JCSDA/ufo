/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSSEAICETHICKNESS_H_
#define UFO_OBSSEAICETHICKNESS_H_

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
#include "util/ObjectCounter.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Total ice concentration observation for UFO.
template <typename MODEL>
class ObsSeaIceThickness : public oops::ObsOperatorBase<MODEL>,
                  private util::ObjectCounter<ObsSeaIceThickness<MODEL>> {
 public:
  static const std::string classname() {return "ufo::ObsSeaIceThickness";}

  ObsSeaIceThickness(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsSeaIceThickness();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperSeaIceThickness_;}
  const int & toFortran() const {return keyOperSeaIceThickness_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperSeaIceThickness_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsSeaIceThickness<MODEL>::ObsSeaIceThickness(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperSeaIceThickness_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_seaicethick_setup_f90(keyOperSeaIceThickness_, &configc);
  const std::vector<std::string> vv{"ice_concentration", "ice_thickness"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsSeaIceThickness created." << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsSeaIceThickness<MODEL>::~ObsSeaIceThickness() {
  ufo_seaicethick_delete_f90(keyOperSeaIceThickness_);
  oops::Log::trace() << "ObsSeaIceThickness destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsSeaIceThickness<MODEL>::obsEquiv(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_seaicethick_eqv_f90(keyOperSeaIceThickness_, gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsSeaIceThickness<MODEL>::print(std::ostream & os) const {
  os << "ObsSeaIceThickness::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSSEAICETHICKNESS_H_
