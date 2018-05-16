/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSAODTLAD_H_
#define UFO_OBSAODTLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperBase.h"
#include "ioda/ObsSpace.h"
#include "util/ObjectCounter.h"
#include "util/Logger.h"

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
/// Aod (currently only temperature) observation for UFO.
template <typename MODEL>
class ObsAodTLAD : public oops::LinearObsOperBase<MODEL>,
                          private util::ObjectCounter<ObsAodTLAD<MODEL>> {
public:
  static const std::string classname() {return "ufo::ObsAodTLAD";}

  ObsAodTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAodTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void obsEquivTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void obsEquivAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

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
ObsAodTLAD<MODEL>::ObsAodTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperAod_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_aod_tlad_setup_f90(keyOperAod_, &configc);
  const std::vector<std::string> vv{"temperature","humidity_mixing_ratio",
      "air_pressure","air_pressure_levels",
      "sulf","bc1","bc2","oc1","oc2","dust1","dust2","dust3","dust4","dust5",
      "seas1","seas2","seas3","seas4","p25"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsAodTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsAodTLAD<MODEL>::~ObsAodTLAD() {
  oops::Log::trace() << "ObsAodTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAodTLAD<MODEL>::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_aod_tlad_settraj_f90(keyOperAod_, geovals.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAodTLAD<MODEL>::obsEquivTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                               const ObsBiasIncrement & bias) const {
  ufo_aod_tlad_eqv_tl_f90(keyOperAod_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAodTLAD<MODEL>::obsEquivAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                               ObsBiasIncrement & bias) const {
  ufo_aod_tlad_eqv_ad_f90(keyOperAod_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsAodTLAD<MODEL>::print(std::ostream & os) const {
  os << "ObsAodTLAD::print not implemented" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSAODTLAD_H_
