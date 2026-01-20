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
#include <memory>
#include <string>
#include "eckit/exception/Exceptions.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathGenerator.h"
#include "ufo/operators/gnssro/utils/StraightLineRayPathGenerator.h"

namespace ufo {

// Factory method for GnssroRayPathGenerator objects.
std::unique_ptr<GnssroRayPathGenerator>
GnssroRayPathGenerator::create(const ioda::ObsSpace & odb,
                               const GnssroRayPathParameters & params)
{
  if (params.rayPathGenType() == StraightLineRayPathGenerator::name()) {
    return std::make_unique<StraightLineRayPathGenerator>(odb, params);
  } else {
    throw eckit::BadValue(
        "Unsupported GnssroRayPathGenerator subclass: " + params.rayPathGenType()
        + "; the only supported type is " + StraightLineRayPathGenerator::name(),
        Here());
  }
}

}  // namespace ufo

