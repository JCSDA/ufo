/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_LAPSERATE_HELPER_H_
#define UFO_PREDICTORS_LAPSERATE_HELPER_H_

#include <map>
#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace ioda {
  class ObsSpace;
}

namespace oops {
  class ObsVariables;
  class Variables;
}

namespace ufo {

class GeoVaLs;
class ObsDiagnostics;

// -----------------------------------------------------------------------------

/// Shared YAML for file or online mean lapse rates (flat under predictor options).
class LapseRateMeanParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LapseRateMeanParameters, Parameters)

 public:
  /// Path to a mean lapse-rate file. Either this must be provided or \ref onlineMean
  /// must be true.
  oops::OptionalParameter<std::string> tlapse{"tlapse", this};

  /// If true, per-channel means are computed online on the first compute for this obs space
  /// (then frozen). Defaults to false. Exactly one of this (when true) and \ref tlapse must be
  /// configured.
  oops::Parameter<bool> onlineMean{"online mean", false, this};

  /// Path to a file to save the per-channel means to.
  oops::OptionalParameter<std::string> saveFile{"save means", this};
};

// -----------------------------------------------------------------------------

class LapseRateHelper {
 public:
  static oops::Variables requiredGeovars();

  LapseRateHelper(const LapseRateMeanParameters & parameters,
                  const oops::ObsVariables & vars);

  oops::ObsVariables requiredHdiags() const;

  /// Fills per-loc tlapchn (row = assim channel index) and shared per-channel means for ObsSpace.
  void fillTlapchnAndMeans(const ioda::ObsSpace & odb,
                           const GeoVaLs & geovals,
                           const ObsDiagnostics & ydiags,
                           std::vector<std::vector<float>> & tlapchn_by_chan,
                           std::vector<float> & means_by_chan) const;

 private:
  void fillTlapchn(const ioda::ObsSpace & odb,
                   const GeoVaLs & geovals,
                   const ObsDiagnostics & ydiags,
                   std::vector<std::vector<float>> & tlapchn_by_chan) const;

  bool onlineMean_;
  std::string tlapseFile_;
  std::string saveFile_;
  std::vector<int> channels_;
  std::map<int, float> fileMeans_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_LAPSERATE_HELPER_H_
