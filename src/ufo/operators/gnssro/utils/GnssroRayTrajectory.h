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
#ifndef UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYTRAJECTORY_H_
#define UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYTRAJECTORY_H_

#include <cassert>
#include <cstddef>  // For std::size_t
#include <map>
#include <vector>

namespace ufo {

// -----------------------------------------------------------------------------
/// Class to represent represent the trajectory of a single GNSSRO ray,
/// used to implement the Tangent Linear and Adjoint models.
/// Each tangent-point sample in a set of RO observations becomes a ray.
class GnssroRayTrajectory
{
 public:
  typedef std::vector<double> DoubleArray_;
  typedef std::vector<int> IntArray_;

  GnssroRayTrajectory();
  GnssroRayTrajectory(const GnssroRayTrajectory & other);
  GnssroRayTrajectory & operator=(const GnssroRayTrajectory & other);
  ~GnssroRayTrajectory();

  // Accessors
  bool isSet() const {return isSet_;}
  std::size_t numNodes() {return wf_.size();}

  // Read-only accessors to array data.
  double jacobianT(uint32_t nodeIdx) const {return jacT_[nodeIdx];}
  double jacobianP(uint32_t nodeIdx) const {return jacP_[nodeIdx];}
  double jacobianQ(uint32_t nodeIdx) const {return jacQ_[nodeIdx];}
  double wf(uint32_t nodeIdx) const {return wf_[nodeIdx];}
  int wi(uint32_t nodeIdx) const {return wi_[nodeIdx];}

  // Read-write accessors to array data.
  double & jacobianT(uint32_t nodeIdx) {assert(nodeIdx < jacT_.size()); return jacT_[nodeIdx];}
  double & jacobianP(uint32_t nodeIdx) {assert(nodeIdx < jacP_.size()); return jacP_[nodeIdx];}
  double & jacobianQ(uint32_t nodeIdx) {assert(nodeIdx < jacQ_.size()); return jacQ_[nodeIdx];}
  double & wf(uint32_t nodeIdx) {assert(nodeIdx < wf_.size()); return wf_[nodeIdx];}
  int & wi(uint32_t nodeIdx) {assert(nodeIdx < wi_.size()); return wi_[nodeIdx];}

  // Methods to call before and after populating the trajectory.
  // Trajectory values should be set using the read-write accessors.
  void initialize(std::size_t numNodes);
  void finalize();

 private:
  //  Member data.
  DoubleArray_ jacT_;  // Derivative of sim var wrt temperature for each node
  DoubleArray_ jacP_;  // Derivative of sim var wrt pressure for each node
  DoubleArray_ jacQ_;  // Derivative of sim var wrt specific humidity for each node
  DoubleArray_ wf_;    // vertical interpolation weight for each node
  IntArray_ wi_;       // vertical interpolation base index (1-based) for each node
  bool isSet_;         // True if trajectory variables that follow are populated.
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_GNSSRORAYTRAJECTORY_H_
