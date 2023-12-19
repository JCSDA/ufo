/*
 * (C) Copyright 2020 NOAA NWS NCEP EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_LAMDOMAINCHECK_LAMDOMAINCHECK_H_
#define UFO_FILTERS_OBSFUNCTIONS_LAMDOMAINCHECK_LAMDOMAINCHECK_H_

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/LAMDomainCheck/LAMDomainCheck.interface.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class LAMDomainCheckParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LAMDomainCheckParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> mapproj{"map_projection", this};
  oops::Parameter<float> esg_a{"a", 0.0f, this};
  oops::Parameter<float> esg_k{"k", 0.0f, this};
  oops::Parameter<float> esg_plat{"plat", 0.0f, this};
  oops::Parameter<float> esg_plon{"plon", 0.0f, this};
  oops::Parameter<float> esg_pazi{"pazi", 0.0f, this};
  oops::Parameter<float> esg_dx{"dx", 1.0f, this};
  oops::Parameter<float> esg_dy{"dy", 1.0f, this};
  oops::Parameter<int> esg_npx{"npx", 2, this};
  oops::Parameter<int> esg_npy{"npy", 2, this};
  oops::Parameter<int> esg_nbdy{"nbdy", 0, this};
  oops::Parameter<bool> save{"save", false, this};

// for a circle domain on sphere with central lat/lon in degree
// and radius in km
  oops::Parameter<float> cenlat{"cenlat", 0.0f, this};
  oops::Parameter<float> cenlon{"cenlon", 0.0f, this};
  oops::Parameter<float> radius{"radius", 100.0f, this};
};

// -----------------------------------------------------------------------------

class LAMDomainCheck : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "LAMDomainCheck";}

  explicit LAMDomainCheck(const eckit::LocalConfiguration &
                                 = eckit::LocalConfiguration());
  ~LAMDomainCheck();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  LAMDomainCheckParameters options_;
};

// -----------------------------------------------------------------------------
}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_LAMDOMAINCHECK_LAMDOMAINCHECK_H_
