/*
 * (C) Copyright
 *
 * Licensed under the terms of the Apache Licence Version 2.0
 * http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_PATHSUM_PATHSUMOPER_H_
#define UFO_OPERATORS_PATHSUM_PATHSUMOPER_H_

#include <array>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variable.h"
#include "oops/base/Variables.h"
#include "ufo/ObsOperatorBase.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {
  class ObsVector;
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
// Parameters for PathSumOper
class PathSumOperParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(PathSumOperParameters, ObsOperatorParametersBase)
/*!
 * \brief Observation operator for path-integrated GeoVaLs
 *
 * The PathSum operator computes the integral (sum) of a GeoVaL variable along a
 * specified path. Two modes are supported:
 *
 *  - **Vertical integration** (`path type: vertical`, default):
 *    For each observation location (latitude, longitude), integrate the chosen
 *    GeoVaL variable vertically with respect to `height_wrt_surface`.
 *    Integration weights can be provided in one of three ways:
 *      - from a GeoVaL variable (`weight variable`),
 *      - from a YAML vector (`weights`),
 *      - computed as vertical distances between adjacent levels (default).
 *    Optionally, integration can be restricted to a user-defined height range
 *    (`height range: [zmin, zmax]`).
 *
 *  - **Slant integration** (`path type: slant`):
 *    Groups observation points by `sequenceNumber` from MetaData. Each group
 *    defines a slant path through 3D space. Integration is performed along each
 *    path, either using GeoVaL-provided weights or by computing segment lengths
 *    from latitude, longitude, and height coordinates.
 *
 * The final integrated value can be scaled by an optional `scaling factor`
 * (default = 1.0), useful for unit conversion (e.g., converting electron density
 * to Total Electron Content).
 *
 * ### YAML example
 *
 * \code{.yaml}
 * obs operator:
 *   name: PathSum
 *   path type: vertical                  # "slant" option is TO BE IMPLEMENTED
 *   geoval variable: electron_density    # required GeoVaL to integrate
 *   weight variable: weightVarName       # optional, from GeoVaLs
 *   weights: [100, 200, 150]             # optional, from YAML vector, only available for vertical
 *                                        # integration for now
 *   height range: [0, 500000]            # optional, height range above surface (unit is determined
 *                                        # by GeoVals)
 *   scaling factor: 1.0e-07              # optional, default = 1.0 
 * \endcode
 *
 * Output:
 *  - The operator writes integrated values into hofx.
 *  - For slant paths, all points sharing the same `sequenceNumber` receive
 *    the same integrated value.
 */

 public:
  // Name of GeoVaL variable to integrate
  oops::RequiredParameter<std::string> geovalVar{"geoval variable", this};

  // Path type: "vertical" or "slant"
  oops::Parameter<std::string> pathType{"path type", "vertical", this};

  // Optional GeoVaL variable for weights
  oops::OptionalParameter<std::string> weightsVar{"weight variable", this};

  // Optional YAML weights
  oops::OptionalParameter<std::vector<double>> weights{"weights", this};

  // Optional height range restriction: [hmin, hmax]
  oops::OptionalParameter<std::vector<double>> heightRange{"height range", this};

  // Optional: whether to interpolate to exact height range boundaries [hmin, hmax].
  //           Only works when heightrange is defined
  oops::Parameter<bool> interpolateBoundaries{"interpolate boundaries", false, this};

  // Optional scaling factor for geoval height unit. True if km is unit of geoval height unit
  // Otherwise, m (SI unit) is assumed
  oops::Parameter<bool> useKmForHeight{"use km for height", false, this};

  // Optional scaling factor for unit conversion or other purposes
  oops::Parameter<float> scalingFactor{"scaling factor", 1.0, this};
};

// -----------------------------------------------------------------------------
// Path types
enum class PathType {VERTICAL, SLANT};

// -----------------------------------------------------------------------------
// PathSum observation operator
class PathSumOper : public ObsOperatorBase {
 public:
  // The type of parameters accepted by the constructor of this operator.
  // This typedef is used by the ObsOperatorFactory.
  typedef PathSumOperParameters Parameters_;

  static const std::string classname() {return "ufo::PathSumOper";}

  PathSumOper(const ioda::ObsSpace &odb, const Parameters_ &params);
  virtual ~PathSumOper();

  // Obs Operator
  void simulateObs(const GeoVaLs &,
                   ioda::ObsVector &,
                   ObsDiagnostics & , const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override {return requiredVars_;}

 private:
  void print(std::ostream &) const override;

  // Compute 3D segment length
  double computeSegmentLength(const std::array<double, 3> &p1,
                              const std::array<double, 3> &p2) const;

  // Compute trapezoidal weight for a given level based on segment lengths
  double computeTrapezoidalWeight(std::size_t lev,
                                  const std::vector<double> &heights,
                                  float lat, float lon,
                                  std::size_t nlevs) const;

  // Interpolate upper/lower boundary if needed
  void interpolateBoundaries(std::vector<double> &heights,
                             std::vector<double> &vals,
                             double hmin, double hmax) const;

  /// Parameters stored
  const ioda::ObsSpace & odb_;
  PathType pathType_;
  const oops::Variable geovalVar_;
  boost::optional<oops::Variable> weightsVar_;
  std::vector<double> weights_;
  std::vector<double> heightRange_;
  bool interpolateBoundaries_;
  bool useKmForHeight_;
  float scalingFactor_;

  oops::Variables requiredVars_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_PATHSUM_PATHSUMOPER_H_
