/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSAIRCRAFTTLAD_H_
#define UFO_OBSAIRCRAFTTLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperBase.h"
#include "ioda/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Logger.h"

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
/// Aircraft (currently only temperature) observation for UFO.
template <typename MODEL>
class ObsAircraftTLAD : public oops::LinearObsOperBase<MODEL>,
                          private util::ObjectCounter<ObsAircraftTLAD<MODEL>> {
public:
  static const std::string classname() {return "ufo::ObsAircraftTLAD";}

  ObsAircraftTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAircraftTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void obsEquivTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void obsEquivAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

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
ObsAircraftTLAD<MODEL>::ObsAircraftTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperAircraft_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_aircraft_tlad_setup_f90(keyOperAircraft_, &configc);
  const std::vector<std::string> vv{"virtual_temperature", "atmosphere_ln_pressure_coordinate"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsAircraftTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsAircraftTLAD<MODEL>::~ObsAircraftTLAD() {
  ufo_aircraft_tlad_delete_f90(keyOperAircraft_);
  oops::Log::trace() << "ObsAircraftTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAircraftTLAD<MODEL>::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_aircraft_tlad_settraj_f90(keyOperAircraft_, geovals.toFortran(), odb_.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAircraftTLAD<MODEL>::obsEquivTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                             const ObsBiasIncrement & bias) const {
  ufo_aircraft_tlad_t_eqv_tl_f90(keyOperAircraft_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAircraftTLAD<MODEL>::obsEquivAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                             ObsBiasIncrement & bias) const {
  ufo_aircraft_tlad_t_eqv_ad_f90(keyOperAircraft_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAircraftTLAD<MODEL>::print(std::ostream & os) const {
  os << "ObsAircraftTLAD::print not implemented" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSAIRCRAFTTLAD_H_
