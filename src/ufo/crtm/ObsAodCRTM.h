/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_CRTM_OBSAODCRTM_H_
#define UFO_CRTM_OBSAODCRTM_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/crtm/ObsAodCRTM.interface.h"
#include "ufo/ObsOperatorBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

// -----------------------------------------------------------------------------
/// AodCRTM observation for UFO.
class ObsAodCRTM : public ObsOperatorBase,
                    private util::ObjectCounter<ObsAodCRTM> {
 public:
  static const std::string classname() {return "ufo::ObsAodCRTM";}

  ObsAodCRTM(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAodCRTM();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}
  const oops::Variables & observed() const {return *varout_;}

  int & toFortran() {return keyOperAodCRTM_;}
  const int & toFortran() const {return keyOperAodCRTM_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperAodCRTM_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
  std::unique_ptr<const oops::Variables> varout_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_CRTM_OBSAODCRTM_H_
