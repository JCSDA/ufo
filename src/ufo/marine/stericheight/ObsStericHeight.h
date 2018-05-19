/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSSTERICHEIGHT_H_
#define UFO_OBSSTERICHEIGHT_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/ObsOperatorBase.h"
#include "ioda/ObsSpace.h"
#include "oops/util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class Locations;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// Steric height/ sea-level observation for UFO.
template <typename MODEL>
class ObsStericHeight : public oops::ObsOperatorBase<MODEL>,
                  private util::ObjectCounter<ObsStericHeight<MODEL>> {
 public:
  static const std::string classname() {return "ufo::ObsStericHeight";}

  ObsStericHeight(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsStericHeight();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperStericHeight_;}
  const int & toFortran() const {return keyOperStericHeight_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperStericHeight_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsStericHeight<MODEL>::ObsStericHeight(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperStericHeight_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_stericheight_setup_f90(keyOperStericHeight_, &configc);
  const std::vector<std::string> vv{"sea_surface_height_above_geoid",
                                    "ocean_potential_temperature",
                                    "ocean_salinity"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsStericHeight created." << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsStericHeight<MODEL>::~ObsStericHeight() {
  ufo_stericheight_delete_f90(keyOperStericHeight_);
  oops::Log::trace() << "ObsStericHeight destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsStericHeight<MODEL>::obsEquiv(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_stericheight_eqv_f90(keyOperStericHeight_, gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsStericHeight<MODEL>::print(std::ostream & os) const {
  os << "ObsStericHeight::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSSTERICHEIGHT_H_
