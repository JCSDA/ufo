/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSINSITUTEMPERATURE_H_
#define UFO_OBSINSITUTEMPERATURE_H_

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
/// Total ice concentration observation for UFO.
template <typename MODEL>
class ObsInsituTemperature : public oops::ObsOperatorBase<MODEL>,
                  private util::ObjectCounter<ObsInsituTemperature<MODEL>> {
 public:
  static const std::string classname() {return "ufo::ObsInsituTemperature";}

  ObsInsituTemperature(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsInsituTemperature();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperInsituTemperature_;}
  const int & toFortran() const {return keyOperInsituTemperature_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperInsituTemperature_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
  ObsInsituTemperature<MODEL>::ObsInsituTemperature(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperInsituTemperature_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_insitutemperature_setup_f90(keyOperInsituTemperature_, &configc);
  const std::vector<std::string> vv{"ocean_potential_temperature", "ocean_salinity", "ocean_layer_thickness"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsInsituTemperature created." << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsInsituTemperature<MODEL>::~ObsInsituTemperature() {
  ufo_insitutemperature_delete_f90(keyOperInsituTemperature_);
  oops::Log::trace() << "ObsInsituTemperature destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
  void ObsInsituTemperature<MODEL>::obsEquiv(const GeoVaLs & gom, ioda::ObsVector & ovec, const ObsBias & bias) const {
  ufo_insitutemperature_eqv_f90(keyOperInsituTemperature_, gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsInsituTemperature<MODEL>::print(std::ostream & os) const {
  os << "ObsInsituTemperature::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSINSITUTEMPERATURE_H_
