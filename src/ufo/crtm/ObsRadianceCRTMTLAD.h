/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_CRTM_OBSRADIANCECRTMTLAD_H_
#define UFO_CRTM_OBSRADIANCECRTMTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/crtm/ObsRadianceCRTMTLAD.interface.h"
#include "ufo/LinearObsOperatorBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
  class LocalConfiguration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// RadianceCRTM (currently only temperature) observation for UFO.
class ObsRadianceCRTMTLAD : public LinearObsOperatorBase,
                        private util::ObjectCounter<ObsRadianceCRTMTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsRadianceCRTMTLAD";}

  ObsRadianceCRTMTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadianceCRTMTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

  // Other
  const oops::Variables & variables() const {return varin_;}

  int & toFortran() {return keyOperRadianceCRTM_;}
  const int & toFortran() const {return keyOperRadianceCRTM_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperRadianceCRTM_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_CRTM_OBSRADIANCECRTMTLAD_H_
