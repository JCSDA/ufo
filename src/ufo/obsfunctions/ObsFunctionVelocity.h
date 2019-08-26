/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSFUNCTIONS_OBSFUNCTIONVELOCITY_H_
#define UFO_OBSFUNCTIONS_OBSFUNCTIONVELOCITY_H_

#include "ufo/obsfunctions/ObsFunctionBase.h"

namespace oops {
  class Variables;
}

namespace ufo {

// -----------------------------------------------------------------------------

class ObsFunctionVelocity : public ObsFunctionBase {
 public:
  ObsFunctionVelocity();
  ~ObsFunctionVelocity();

  void compute(const ioda::ObsDataVector<float> &,
               const ioda::ObsDataVector<float> &,
               ioda::ObsDataVector<float> &) const;
  const oops::Variables & requiredObsData() const;
  const oops::Variables & requiredMetaData() const;
 private:
  oops::Variables obsvars_;
  oops::Variables metadatavars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSFUNCTIONS_OBSFUNCTIONVELOCITY_H_
