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
#ifndef UFO_OPERATORS_GNSSRO_UTILS_GNSSROPROFILEEXTRACTOR_H_
#define UFO_OPERATORS_GNSSRO_UTILS_GNSSROPROFILEEXTRACTOR_H_

#include <cstddef>  // For std::size_t
#include <vector>
#include "ioda/ObsSpace.h"
#include "ufo/operators/gnssro/utils/GnssroProfileSlice.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Class to extract the GnssroProfileSlices from a set of GNSSRO observation samples.
/// Each slice references a distinct RO profile.
class GnssroProfileExtractor {
 public:
  typedef std::vector<GnssroProfileSlice> Slices;
  typedef Slices::const_iterator const_iterator;

  explicit GnssroProfileExtractor(const ioda::ObsSpace & odb);
  ~GnssroProfileExtractor();

  const Slices & slices() const {return slices_;}
  const_iterator cbegin() const {return slices_.cbegin();}
  const_iterator cend() const {return slices_.cend();}
  bool is_empty() const {return slices_.empty();}

 private:
  GnssroProfileExtractor() = delete;
  GnssroProfileExtractor(const GnssroProfileExtractor &) = delete;
  GnssroProfileExtractor(GnssroProfileExtractor &&) noexcept = delete;
  GnssroProfileExtractor& operator=(const GnssroProfileExtractor &) = delete;
  GnssroProfileExtractor& operator=(GnssroProfileExtractor &&) noexcept = delete;

  Slices slices_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_GNSSROPROFILEEXTRACTOR_H_
