/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSFUNCTIONS_OBSFUNCTION_H_
#define UFO_OBSFUNCTIONS_OBSFUNCTION_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "ioda/ObsDataVector.h"
#include "ufo/obsfunctions/ObsFunctionBase.h"

namespace oops {
  class Variables;
}

namespace ufo {

class GeoVaLs;

// -----------------------------------------------------------------------------

class ObsFunction : private boost::noncopyable {
 public:
/// constructor takes function name (for factory) on input
  explicit ObsFunction(const std::string &);
  ~ObsFunction();

/// compute(metadata, obs values, output)
  void compute(const ioda::ObsDataVector<float> &,
               const ioda::ObsDataVector<float> &,
               const GeoVaLs &,
               ioda::ObsDataVector<float> &) const;
/// required variables (@ObsValue/@HofX)
  const oops::Variables & requiredObsData() const;
/// required metadata
  const oops::Variables & requiredMetaData() const;
/// required geovals
  const oops::Variables & requiredGeoVaLs() const;
 private:
  std::unique_ptr<ObsFunctionBase> obsfct_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSFUNCTIONS_OBSFUNCTION_H_
