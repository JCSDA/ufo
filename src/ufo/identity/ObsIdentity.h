/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_IDENTITY_OBSIDENTITY_H_
#define UFO_IDENTITY_OBSIDENTITY_H_

#include <ostream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/identity/ObsIdentity.interface.h"
#include "ufo/ObsOperatorBase.h"

/// Forward declarations
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
/// Generic identity observation operator class
class ObsIdentity : public ObsOperatorBase,
                          private util::ObjectCounter<ObsIdentity> {
 public:
  static const std::string classname() {return "ufo::ObsIdentity";}

  ObsIdentity(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsIdentity();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}
  const oops::Variables & observed() const {return *varout_;}

  int & toFortran() {return keyOperObsIdentity_;}
  const int & toFortran() const {return keyOperObsIdentity_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperObsIdentity_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
  boost::scoped_ptr<const oops::Variables> varout_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_IDENTITY_OBSIDENTITY_H_
