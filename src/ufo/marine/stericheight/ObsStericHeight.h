/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_MARINE_STERICHEIGHT_OBSSTERICHEIGHT_H_
#define UFO_MARINE_STERICHEIGHT_OBSSTERICHEIGHT_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/marine/stericheight/ObsStericHeight.interface.h"
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
  class Locations;
  class ObsBias;

// -----------------------------------------------------------------------------
/// Steric height/ sea-level observation operator class
class ObsStericHeight : public ObsOperatorBase,
                        private util::ObjectCounter<ObsStericHeight> {
 public:
  static const std::string classname() {return "ufo::ObsStericHeight";}

  ObsStericHeight(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsStericHeight();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}
  const oops::Variables & observed() const {return *varout_;}
  Locations * locateObs(const util::DateTime &, const util::DateTime &) const;

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
  boost::scoped_ptr<const oops::Variables> varout_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_MARINE_STERICHEIGHT_OBSSTERICHEIGHT_H_
