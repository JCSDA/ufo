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
#ifndef UFO_OPERATORS_GNSSRO_UTILS_GNSSROPROFILESLICE_H_
#define UFO_OPERATORS_GNSSRO_UTILS_GNSSROPROFILESLICE_H_

#include <cstddef>  // For std::size_t
#include "oops/util/DateTime.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Class to represent represent the boundaries in the samples of observed
/// GNSSRO data between separate, largely vertical profiles.
/// These slices are created by the GnssroProfileExtractor class.
///
class GnssroProfileSlice {
 public:
  GnssroProfileSlice()
    : seqNum_(-1)
    , time_(19000101, 0)
    , start_(0)
    , end_(0) {}
  GnssroProfileSlice(int seqNum, const util::DateTime& time, size_t start, size_t end)
    : seqNum_(seqNum)
    , time_(time)
    , start_(start)
    , end_(end) {}
  GnssroProfileSlice(const GnssroProfileSlice & other)
    : seqNum_(other.seqNum_)
    , time_(other.time_)
    , start_(other.start_)
    , end_(other.end_) {}
  ~GnssroProfileSlice() {}
  GnssroProfileSlice & operator=(const GnssroProfileSlice & other) {
    if (this != &other) {
      seqNum_ = other.seqNum_;
      time_   = other.time_;
      start_  = other.start_;
      end_    = other.end_;
    }
    return *this;
  }

  //  Accessors.
  int seqNum() const {return seqNum_;}
  std::size_t start() const {return start_;}
  util::DateTime time() const {return time_;}
  std::size_t end() const {return end_;}
  std::size_t size() const {return end_ - start_;}

 private:
  // A slice includes all indices start_ <= idx < end_,
  // i.e. inclusive of start and exclusive of end.
  int            seqNum_;
  util::DateTime time_;
  std::size_t    start_;
  std::size_t    end_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_GNSSROPROFILESLICE_H_
