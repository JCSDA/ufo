/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSAIRCRAFT_H_
#define UFO_OBSAIRCRAFT_H_

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

namespace ufo {

// -----------------------------------------------------------------------------
/// Aircraft (currently only temperature) observation for UFO.
template <typename MODEL>
class ObsAircraft : public oops::ObsOperatorBase<MODEL>,
                  private util::ObjectCounter<ObsAircraft<MODEL>> {
 public:
  static const std::string classname() {return "ufo::ObsAircraft";}

  ObsAircraft(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAircraft();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperAircraft_;}
  const int & toFortran() const {return keyOperAircraft_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperAircraft_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsAircraft<MODEL>::ObsAircraft(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperAircraft_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_aircraft_setup_f90(keyOperAircraft_, &configc);
  const std::vector<std::string> vv{"virtual_temperature", "atmosphere_ln_pressure_coordinate"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsAircraft created." << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsAircraft<MODEL>::~ObsAircraft() {
  ufo_aircraft_delete_f90(keyOperAircraft_);
  oops::Log::trace() << "ObsAircraft destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAircraft<MODEL>::obsEquiv(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_aircraft_t_eqv_f90(keyOperAircraft_, gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAircraft<MODEL>::print(std::ostream & os) const {
  os << "ObsAircraft::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSAIRCRAFT_H_
