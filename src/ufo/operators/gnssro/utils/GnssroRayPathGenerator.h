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
#ifndef UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATHGENERATOR_H_
#define UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATHGENERATOR_H_

#include <memory>
#include <string>

#include "ioda/ObsSpace.h"
#include "ufo/operators/gnssro/utils/GnssroProfileRayPaths.h"
#include "ufo/operators/gnssro/utils/GnssroProfileSlice.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathParameters.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Base class for GnssroRayPathGenerator
///
/// The developer should define a subclass of GnssroRayPathGenerator that
/// implements the makeRayPath method.
///
class GnssroRayPathGenerator {
 public:
  static std::string name() {return "GnssroRayPathGenerator";}

  //  Factory method.
  static std::unique_ptr<GnssroRayPathGenerator> create(
      const ioda::ObsSpace & odb, const GnssroRayPathParameters & params);

  //  Constructor/Destructor
  explicit GnssroRayPathGenerator(const ioda::ObsSpace & odb,
                                  const GnssroRayPathParameters & params)
     : odb_(odb)
     , params_(params) {}
  virtual ~GnssroRayPathGenerator() {}

  /// \brief Return a GnssroRayPath object based on observation data read for the specified
  /// GNSSRO profile from the ObsSpace.
  virtual std::unique_ptr<GnssroProfileRayPaths> makeProfileRayPaths(
      const GnssroProfileSlice & profileSlice) = 0;

 protected:
  const ioda::ObsSpace& odb_;
  GnssroRayPathParameters params_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYPATHGENERATOR_H_
