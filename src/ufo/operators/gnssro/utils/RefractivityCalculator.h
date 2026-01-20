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
#ifndef UFO_OPERATORS_GNSSRO_UTILS_REFRACTIVITYCALCULATOR_H_
#define UFO_OPERATORS_GNSSRO_UTILS_REFRACTIVITYCALCULATOR_H_

#include <memory>
#include <string>

namespace ufo {

// -----------------------------------------------------------------------------
/// Base class for RefractivityCalculator.
///
/// Refractivity algorithms are identified by a textual name. Each algorithm
/// is implemented as a separate subclass of RefractivityCalculator.
/// In addition, the same algorithm can be identified as using a compressible
/// or non-compressible atmosphere using the useCompress argument.
///
/// The RefractivityCalculator can compute refractivity from air temperature
/// (T) in Kelvin, specific humidity (Q) in kg/kg and pressure (P) in Pa.
/// It can also compute the derivatives of refractivity with respect to each
/// of these variables.
///
class RefractivityCalculator {
 public:
  //  Refractivity units support:
  //  Conversions between excess index of refraction (n-1) to Refractivity N.
  static constexpr double EXCESS_IOR_TO_N = 1e6;
  static constexpr double N_TO_EXCESS_IOR = 1e-6;

  //  Factory method.
  static std::unique_ptr<RefractivityCalculator> create(
      const std::string& algorithmName, bool useCompress);

  //  Constructor/Destructor
  RefractivityCalculator() = default;
  virtual ~RefractivityCalculator() = default;
  virtual const char * algorithmName() const = 0;

  /// \brief Return refractivity as a function of temperature, specific
  /// humidity, and pressure.
  virtual double N(double T, double Q, double P) const = 0;

  /// \brief Return derivative of refractivity with respect to temperature
  virtual double dNdT(double T, double Q, double P) const = 0;
  /// \brief Return derivative of refractivity with respect to specific humidity
  virtual double dNdQ(double T, double Q, double P) const = 0;
  /// \brief Return derivative of refractivity with respect to pressure
  virtual double dNdP(double T, double Q, double P) const = 0;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_UTILS_REFRACTIVITYCALCULATOR_H_
