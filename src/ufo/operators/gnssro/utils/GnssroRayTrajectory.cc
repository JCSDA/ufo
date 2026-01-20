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
#include <algorithm>
#include <limits>
#include <numeric>
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/operators/gnssro/utils/GnssroRayTrajectory.h"

namespace ufo {

// -----------------------------------------------------------------------------
GnssroRayTrajectory::GnssroRayTrajectory()
  : jacT_()
  , jacP_()
  , jacQ_()
  , wf_()
  , wi_()
  , isSet_(false)
{
}

// -----------------------------------------------------------------------------
GnssroRayTrajectory::GnssroRayTrajectory(const GnssroRayTrajectory& other)
  : jacT_(other.jacT_)
  , jacP_(other.jacP_)
  , jacQ_(other.jacQ_)
  , wf_(other.wf_)
  , wi_(other.wi_)
  , isSet_(other.isSet_)
{
}

// -----------------------------------------------------------------------------
GnssroRayTrajectory & GnssroRayTrajectory::operator=(const GnssroRayTrajectory& other)
{
  if (this != &other)
  {
    jacT_ = other.jacT_;
    jacP_ = other.jacP_;
    jacQ_ = other.jacQ_;
    wf_ = other.wf_;
    wi_ = other.wi_;
    isSet_ = other.isSet_;
  }
  return *this;
}

// -----------------------------------------------------------------------------
GnssroRayTrajectory::~GnssroRayTrajectory()
{
}

// -----------------------------------------------------------------------------
void GnssroRayTrajectory::initialize(std::size_t numNodes)
{
  isSet_ = false;
  jacT_.assign(numNodes, util::missingValue<double>());
  jacP_.assign(numNodes, util::missingValue<double>());
  jacQ_.assign(numNodes, util::missingValue<double>());
  wf_.assign(numNodes, util::missingValue<double>());
  wi_.assign(numNodes, 0);
}

// -----------------------------------------------------------------------------
void GnssroRayTrajectory::finalize()
{
  isSet_ = true;
  return;
}
// -----------------------------------------------------------------------------

}  // namespace ufo

