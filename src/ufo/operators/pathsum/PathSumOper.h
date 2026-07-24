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
 * The pathsum operator computes a weighted summation of a GeoVaLs variable
 * along a vertical or slant path. When the weights correspond to path lengths,
 * the summation provides a numerical approximation to the path integral.

 * Integration weights may be:
 *   - read from GeoVaLs
 *   - specified directly in YAML
 *   - computed internally using trapezoidal integration
 *
 * The operator also supports:
 *   - optional height-range restriction
 *   - optional interpolation to exact height boundaries at two the end points of
 *     the path, implemented only for vertical paths
 *
 * ### YAML example
 *
 * \code{.yaml}
 * obs operator:
 *   name: PathSum
 *   path type: vertical                  
 *   geoval variable: electron_density    # GeoVaLs variable to integrate
 *   weights: [100, 200, 150]             # optional, weights from YAML
 *   height range: [0, 500000]            # optional, height range 
 *                                        # (height unit is determined by GeoVaLs)
 *   scaling factor: 1.0e-07              # optional, default = 1.0
 * \endcode
 *
 * Output:
 *  - The operator writes integrated values into HofX.
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

  // Optional GeoVaL variables for slant path coordinates
  oops::Parameter<std::string> pathPointLatVar{"path point latitude variable",
                                               "pathPointLatitude", this};
  oops::Parameter<std::string> pathPointLonVar{"path point longitude variable",
                                               "pathPointLongitude", this};
  // Optional GeoVaL variables for vertical or slant path coordinates
  oops::Parameter<std::string> pathPointHeightVar{"path point height variable",
                                                  "pathPointHeight", this};
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

  // Geometry container for a path
  struct PathGeometry {
    std::vector<float> lat;
    std::vector<float> lon;
    std::vector<double> height;
  };

  // Compute 3D segment length
  double computeSegmentLength(const std::array<double, 3> &p1,
                              const std::array<double, 3> &p2) const;

  // Compute trapezoidal weight for a given ipoint based on segment lengths
  double computeTrapezoidalWeight(std::size_t ipoint,
                                  const PathGeometry &geom) const;

  // Interpolate upper/lower boundary if needed
  void interpolateBoundaries(std::vector<double> &heights,
                             std::vector<double> &vals,
                             double hmin, double hmax) const;

    // Build geometry for vertical or slant paths
  PathGeometry buildGeometry(std::size_t loc,
                             const GeoVaLs &geovals,
                             const std::vector<float> &lats,
                             const std::vector<float> &lons,
                             bool needHeights,
                             bool needPathLatLon) const;

  // Build weights
  void buildWeights(std::vector<double> &wts,
                    const PathGeometry &geom,
                    std::size_t loc,
                    const GeoVaLs &geovals,
                    bool useHeightRange,
                    const double hmin, const double hmax,
                    const std::size_t npoint) const;

  // Generic integration
  double integrate(const std::vector<double> &vals,
                   const std::vector<double> &wts) const;

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
  boost::optional<oops::Variable> pathPointLatVar_;
  boost::optional<oops::Variable> pathPointLonVar_;
  boost::optional<oops::Variable> pathPointHeightVar_;

  oops::Variables requiredVars_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_PATHSUM_PATHSUMOPER_H_
