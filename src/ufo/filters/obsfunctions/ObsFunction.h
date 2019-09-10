/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTION_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTION_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"


namespace ufo {

// -----------------------------------------------------------------------------

class ObsFunction : private boost::noncopyable {
 public:
/// constructor takes function name (for factory) on input
  explicit ObsFunction(const std::string &);
  ~ObsFunction();

/// compute(metadata, obs values, output)
  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
/// required variables
  const ufo::Variables & requiredVariables() const;
 private:
  std::unique_ptr<ObsFunctionBase> obsfct_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTION_H_
