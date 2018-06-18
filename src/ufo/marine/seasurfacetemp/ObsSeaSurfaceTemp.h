/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSSEASURFACETEMP_H_
#define UFO_OBSSEASURFACETEMP_H_

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
class ObsSeaSurfaceTemp : public oops::ObsOperatorBase<MODEL>,
                  private util::ObjectCounter<ObsSeaSurfaceTemp<MODEL>> {
 public:
  static const std::string classname() {return "ufo::ObsSeaSurfaceTemp";}

  ObsSeaSurfaceTemp(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsSeaSurfaceTemp();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperSeaSurfaceTemp_;}
  const int & toFortran() const {return keyOperSeaSurfaceTemp_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperSeaSurfaceTemp_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
  ObsSeaSurfaceTemp<MODEL>::ObsSeaSurfaceTemp(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperSeaSurfaceTemp_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_seasurfacetemp_setup_f90(keyOperSeaSurfaceTemp_, &configc);
  const std::vector<std::string> vv{"ocean_upper_level_temperature"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsSeaSurfaceTemp created." << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsSeaSurfaceTemp<MODEL>::~ObsSeaSurfaceTemp() {
  ufo_seasurfacetemp_delete_f90(keyOperSeaSurfaceTemp_);
  oops::Log::trace() << "ObsSeaSurfaceTemp destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
  void ObsSeaSurfaceTemp<MODEL>::obsEquiv(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_seasurfacetemp_eqv_f90(keyOperSeaSurfaceTemp_, gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsSeaSurfaceTemp<MODEL>::print(std::ostream & os) const {
  os << "ObsSeaSurfaceTemp::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSSEASURFACETEMP_H_
