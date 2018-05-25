/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSINSITUTEMPERATURETLAD_H_
#define UFO_OBSINSITUTEMPERATURETLAD_H_

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
/// Temperature Profile observation for model.
template <typename MODEL>
class ObsInsituTemperatureTLAD : public oops::LinearObsOperBase<MODEL>, 
                              private util::ObjectCounter<ObsInsituTemperatureTLAD<MODEL>> {
public:
  static const std::string classname() {return "ufo::ObsInsituTemperatureTLAD";}

  ObsInsituTemperatureTLAD(const ioda::ObsSpace &, const eckit::Configuration &);    
  virtual ~ObsInsituTemperatureTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void obsEquivTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void obsEquivAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

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
  ObsInsituTemperatureTLAD<MODEL>::ObsInsituTemperatureTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperInsituTemperature_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_insitutemperature_tlad_setup_f90(keyOperInsituTemperature_, &configc);
  const std::vector<std::string> vv{"ocean_potential_temperature", "ocean_salinity", "ocean_layer_thickness"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsInsituTemperatureTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsInsituTemperatureTLAD<MODEL>::~ObsInsituTemperatureTLAD() {
  ufo_insitutemperature_tlad_delete_f90(keyOperInsituTemperature_);
  oops::Log::trace() << "ObsInsituTemperatureTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsInsituTemperatureTLAD<MODEL>::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_insitutemperature_tlad_settraj_f90(keyOperInsituTemperature_, geovals.toFortran(), odb_.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
  void ObsInsituTemperatureTLAD<MODEL>::obsEquivTL(const GeoVaLs & geovals, ioda::ObsVector & ovec, const ObsBiasIncrement & bias) const {
  ufo_insitutemperature_tlad_eqv_tl_f90(keyOperInsituTemperature_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
  void ObsInsituTemperatureTLAD<MODEL>::obsEquivAD(GeoVaLs & geovals, const ioda::ObsVector & ovec, ObsBiasIncrement & bias) const {
  ufo_insitutemperature_tlad_eqv_ad_f90(keyOperInsituTemperature_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsInsituTemperatureTLAD<MODEL>::print(std::ostream & os) const {
  os << "ObsInsituTemperatureTLAD::print not implemented" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSINSITUTEMPERATURETLAD_H_
