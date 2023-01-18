/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OCEANPRESSURETODEPTH_H_
#define UFO_FILTERS_OBSFUNCTIONS_OCEANPRESSURETODEPTH_H_

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ObsFilterData;

/// \brief Options controlling OceanPressureToDepth ObsFunction
class OceanPressureToDepthParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OceanPressureToDepthParameters, Parameters)

 public:
  /// Input variables of the ocean pressure to depth conversion
  oops::RequiredParameter<Variable> pressure{"pressure variable", this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs depth converted from pressure (reference not given in OPS, but
///        appears to be from https://archimer.ifremer.fr/doc/00447/55889/57949.pdf )
///
/// Example
///
///  obs function:
///    name: ObsFunction/OceanPressureToDepth
///    options:
///      pressure variable: ObsValue/waterPressure
///
/// will return depth (m), given pressure (Pa) and latitude (deg). Depth calculated
/// according to formula (3) in https://archimer.ifremer.fr/doc/00447/55889/57949.pdf
///

class OceanPressureToDepth : public ObsFunctionBase<float> {
 public:
  explicit OceanPressureToDepth(const eckit::LocalConfiguration &);
  ~OceanPressureToDepth();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  OceanPressureToDepthParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OCEANPRESSURETODEPTH_H_
