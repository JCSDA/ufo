/*
 * (C) Copyright 2025 Space Sciences and Engineering, LLC (dba PlanetiQ).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * author Steve Marshall (smarshall@planetiq.com)
 */
#ifndef UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATHORCHESTRATOR_H_
#define UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATHORCHESTRATOR_H_

#include <memory>
#include <string>
#include <vector>
#include "ioda/ObsSpace.h"
#include "oops/util/DateTime.h"
#include "oops/util/Range.h"
#include "ufo/operators/gnssro/utils/GnssroProfileExtractor.h"
#include "ufo/operators/gnssro/utils/GnssroProfileRayPaths.h"
#include "ufo/operators/gnssro/utils/GnssroProfileSlice.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathGenerator.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathParameters.h"

namespace ufo {
// Forward class references
class RefractivityCalculator;
class GnssroGeoVaLs;

// -----------------------------------------------------------------------------
/// GnssroRayPathOrchestrator
///
class GnssroRayPathOrchestrator {
 public:
  // Type definitions
  typedef std::unique_ptr<GnssroRayPathGenerator> Generator_;
  typedef std::unique_ptr<GnssroProfileRayPaths> Profile_;
  typedef std::vector<Profile_> Profiles_;
  typedef GnssroProfileRayPaths::Ray_ Ray_;

  // Nested class for holding elements of a trajectory for one
  // node of a 2D ray path.
  class TrajTuple
  {
   public:
    TrajTuple()
      : dNdT_(0.0)
      , dNdQ_(0.0)
      , dNdP_(0.0)
      , wf_(0.0)
      , wi_(0)
    {}

    TrajTuple(double dNdT, double dNdQ, double dNdP, double wf, int wi)
      : dNdT_(dNdT)
      , dNdQ_(dNdQ)
      , dNdP_(dNdP)
      , wf_(wf)
      , wi_(wi)
    {}

    // Read-only accessors.
    double dNdT() const {return dNdT_;}
    double dNdQ() const {return dNdQ_;}
    double dNdP() const {return dNdP_;}
    double wf() const {return wf_;}
    int wi() const {return wi_;}

   private:
    // Member data.
    double dNdT_;
    double dNdQ_;
    double dNdP_;
    double wf_;
    int    wi_;
  };  // end of nested class definition for TrajTuple

  /// \brief Class method to apply pressure interpolation.
  /// This is consistent with interpolation weights and indices
  /// produced by vert_interp_f90. Returns interpolated value
  /// from pres using base index wi (1-based) and weight wf
  /// as well as information from the temperature profile.
  static double pressureInterpApply(
      int nlevs,            // Number of vertical levels
      const double * pres,  // Column of pressure values
      const double * temp,  // Column of temperature values
      int wi,               // Vertical interp base index (1-based)
      double wf);           // Linear vertical interp factor

  /// \brief Class method for computing model refractivity.
  static double modelRefractivity(
      float geomHgt, float lat, const GnssroGeoVaLs & ggv, std::size_t profileIdx,
      const std::unique_ptr<ufo::RefractivityCalculator> & refrCalc);

  /// \brief Class method for computing model trajectory.
  static TrajTuple modelTrajectory(
      float geomHgt, float lat, const GnssroGeoVaLs & ggv, std::size_t profileIdx,
      const std::unique_ptr<ufo::RefractivityCalculator> & refrCalc);

  //  Constructor/Destructor
  explicit GnssroRayPathOrchestrator(
      const ioda::ObsSpace & odb,
      const GnssroRayPathParameters & params);
  ~GnssroRayPathOrchestrator();
  GnssroRayPathOrchestrator() = delete;
  GnssroRayPathOrchestrator(const GnssroRayPathOrchestrator &) = delete;
  GnssroRayPathOrchestrator(GnssroRayPathOrchestrator &&) noexcept = delete;
  GnssroRayPathOrchestrator& operator=(const GnssroRayPathOrchestrator &) = delete;
  GnssroRayPathOrchestrator& operator=(GnssroRayPathOrchestrator &&) noexcept = delete;

  //  Read-only Accessors
  const GnssroRayPathParameters & parameters() const {return params_;}
  const GnssroProfileExtractor & extractor() const {return extractor_;}
  const Generator_ & generator() const {return generator_;}
  const Profiles_ & profiles() const {return profiles_;}

  //  Number of paths (vertical columns) requested from the model.
  std::size_t numSampledLocs() const;

  // Method to allocate and fill the location and time arrays for the SampledLocations.
  // It clears and reallocates the vectors to the necessary number of sampled locations.
  // It returns the common allocation size used for all three vectors.
  std::size_t getSampledLatLonTimes(
          std::vector<float>& lats,
          std::vector<float>& lons,
          std::vector<util::DateTime>& times) const;

  /// Sets the range of SampledLocations used by each ray path.
  void fillPathsGroupedByLocation(
          std::size_t numSampledLocs,
          std::vector<util::Range<std::size_t>> & pathsGroupedByLocation) const;

  //  Number of rays computed within profile ray paths.
  std::size_t totalNumRays() const;

  //  Read-write Accessors
  Generator_ & generator() {return generator_;}
  Profiles_ & profiles() {return profiles_;}

  // ostream operator can access private members.
  friend std::ostream& operator<<(std::ostream& os, const GnssroRayPathOrchestrator& orc);

 private:
  //  Member data
  const ioda::ObsSpace& odb_;
  GnssroRayPathParameters params_;
  GnssroProfileExtractor extractor_;
  Generator_ generator_;
  Profiles_ profiles_;
  std::vector<float> sampledLats_;
  std::vector<float> sampledLons_;
  std::vector<util::DateTime> sampledTimes_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATHORCHESTRATOR_H_
