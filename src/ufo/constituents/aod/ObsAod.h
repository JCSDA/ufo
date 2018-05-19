/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSAOD_H_
#define UFO_OBSAOD_H_

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
/// Aod observation for UFO.
template <typename MODEL>
class ObsAod : public oops::ObsOperatorBase<MODEL>,
                  private util::ObjectCounter<ObsAod<MODEL>> {
 public:
  static const std::string classname() {return "ufo::ObsAod";}

  ObsAod(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAod();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperAod_;}
  const int & toFortran() const {return keyOperAod_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperAod_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsAod<MODEL>::ObsAod(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperAod_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_aod_setup_f90(keyOperAod_, &configc);
  const std::vector<std::string> vv{"temperature","humidity_mixing_ratio",
      "air_pressure","air_pressure_levels",
      "sulf","bc1","bc2","oc1","oc2","dust1","dust2","dust3","dust4","dust5",
      "seas1","seas2","seas3","seas4","p25"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsAod created." << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsAod<MODEL>::~ObsAod() {
  ufo_aod_delete_f90(keyOperAod_);
  oops::Log::trace() << "ObsAod destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAod<MODEL>::obsEquiv(const GeoVaLs & gom, ioda::ObsVector & ovec,
                         const ObsBias & bias) const {
  ufo_aod_eqv_f90(keyOperAod_, gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAod<MODEL>::print(std::ostream & os) const {
  os << "ObsAod::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSAOD_H_
