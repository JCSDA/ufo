/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

// TODO(anyone): through the file replace UFO_EXAMPLE_OBSEXAMPLE_H with the unique string
// (e.g. UFO_EXAMPLE_<YOUR_OBS_OPERATOR_NAME>_H
#ifndef UFO_EXAMPLE_OBSEXAMPLE_H_
#define UFO_EXAMPLE_OBSEXAMPLE_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"

/// Forward declarations
namespace eckit {
class Configuration;
}

namespace ioda {
class ObsVector;
}

namespace ufo {
class GeoVaLs;
class ObsBias;

// -----------------------------------------------------------------------------
/// Example for the observation operator class.
// TODO(anyone): through the file replace ObsExample with <Your_Obs_Operator_Name>
class ObsExample : public ObsOperatorBase,
                   private util::ObjectCounter<ObsExample> {
 public:
  static const std::string classname() {return "ufo::ObsExample";}

  ObsExample(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsExample();

  // Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_EXAMPLE_OBSEXAMPLE_H_
