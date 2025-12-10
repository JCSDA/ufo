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
#include <iomanip>
#include <iostream>
#include <memory>
#include "eckit/exception/Exceptions.h"
#include "ufo/operators/gnssro/utils/RefractivityCalculator.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------
// Specific refractivity algorithms are all declared and implemented in this
// file. The only intended interface if through the RefractivityCalculator
// base class, constructed through the create() factory method.
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Rueger, 2002 (compressible) or Bevis, 1994 (non-compressible)
// -----------------------------------------------------------------------------
class RuegerBevisRefractivityCalculator : public RefractivityCalculator
{
 public:
  static constexpr const char ALGORITHM_NAME[] = "RuegerBevis";
  explicit RuegerBevisRefractivityCalculator(bool useCompress = true);
  ~RuegerBevisRefractivityCalculator() override;
  const char * algorithmName() const override;
  double N(double T, double Q, double P) const override;
  double dNdT(double T, double Q, double P) const override;
  double dNdQ(double T, double Q, double P) const override;
  double dNdP(double T, double Q, double P) const override;
 private:
  double n_a_;
  double n_b_;
  double n_c_;
  double oneMinusRdOverRv_;
};
// Implementation of RuegerBevisRefractivityCalculator

RuegerBevisRefractivityCalculator::RuegerBevisRefractivityCalculator(bool useCompress)
  : RefractivityCalculator()
  , n_a_(0.776890)   // Default compressible, Rueger 2002
  , n_b_(3.75463e3)
  , n_c_(0.712952)
  , oneMinusRdOverRv_(1.0 - Constants::rd_over_rv)
{
  // This combination of refractivity algorithm used in cucurull 2010, Healy 2011
  // Default is for compressible atmosphere, per Rueger 2002.
  if (!useCompress) {
    // Non-compressible atmosphere, per Bevis et al 1994
    n_a_ = 0.7760;
    n_b_ = 3.739e3;
    n_c_ = 0.704;
  }
  n_c_ -= n_a_;
}

RuegerBevisRefractivityCalculator::~RuegerBevisRefractivityCalculator()
{
}

const char * RuegerBevisRefractivityCalculator::algorithmName() const
{
    return ALGORITHM_NAME;
}

double RuegerBevisRefractivityCalculator::N(double T, double Q, double P) const
{
  double tfact = oneMinusRdOverRv_ * Q + Constants::rd_over_rv;
  double refr1 = n_a_ * P / T;
  double refr2 = n_b_ * Q * P / (T * T * tfact);
  double refr3 = n_c_ * Q * P / (T * tfact);
  double refr  = refr1 + refr2 + refr3;
  return refr;
}

double RuegerBevisRefractivityCalculator::dNdT(double T, double Q, double P) const
{
  // Algorithm taken from Jacobian calculation in ufo_gnssro_refncep_tlad_mod.F90
  // self%jac_t(iobs) = - n_a*gesP/gesT**2                                                &
  //                    - n_b*two*gesP*gesQ/( ((1-rd_over_rv)*gesQ+rd_over_rv)*gesT**3  ) &
  //                    - n_c*gesP*gesQ/( ((1-rd_over_rv)*gesQ+rd_over_rv)*gesT**2 )
  double T_sq = T * T;
  double PxQ = P * Q;
  double tfact = oneMinusRdOverRv_ * Q + Constants::rd_over_rv;
  double jacT = - n_a_ * P / T_sq
                - n_b_ * 2.0 * PxQ / (tfact * T_sq * T)
                - n_c_ * PxQ / (tfact * T_sq);
  return jacT;
}

double RuegerBevisRefractivityCalculator::dNdQ(double T, double Q, double P) const
{
  // Algorithm taken from Jacobian calculation in ufo_gnssro_refncep_tlad_mod.F90
  // self%jac_q(iobs)   =   n_b*gesP/( gesT**2*( (1-rd_over_rv)*gesQ+rd_over_rv)**2 )         &
  //                                * rd_over_rv                                              &
  //                      + n_c*gesP/( gesT   *( (1-rd_over_rv)*gesQ+rd_over_rv)**2 )         &
  //                                * rd_over_rv
  double T_sq = T * T;
  double tfact = oneMinusRdOverRv_ * Q + Constants::rd_over_rv;
  double tfact_sq = tfact * tfact;
  double jacQ = n_b_ * P / (T_sq * tfact_sq) * Constants::rd_over_rv
              + n_c_ * P / (T * tfact_sq) * Constants::rd_over_rv;
  return jacQ;
}

double RuegerBevisRefractivityCalculator::dNdP(double T, double Q, double P) const
{
  // Algorithm taken from Jacobian calculation in ufo_gnssro_refncep_tlad_mod.F90
  // self%jac_prs(iobs) = n_a/gesT                                                 &
  //                    + n_b*gesQ/ ( ((1-rd_over_rv)*gesQ+rd_over_rv)*gesT**2 )   &
  //                    + n_c*gesQ/ ( ((1-rd_over_rv)*gesQ+rd_over_rv)*gesT )
  double tfact = oneMinusRdOverRv_ * Q + Constants::rd_over_rv;
  double jacP = n_a_ / T
              + n_b_ * Q / (tfact * T * T)
              + n_c_ * Q / (tfact * T);
  return jacP;
}

// -----------------------------------------------------------------------------
// Base class implementation.
// -----------------------------------------------------------------------------
// Factory method.
std::unique_ptr<RefractivityCalculator> RefractivityCalculator::create(
        const std::string& algorithmName, bool useCompress)
{
  if (algorithmName == RuegerBevisRefractivityCalculator::ALGORITHM_NAME) {
    return std::make_unique<RuegerBevisRefractivityCalculator>(useCompress);
  } else {
    throw eckit::BadValue(
        "Unsupported algorithm name for RefractivityCalculator: " + algorithmName,
        Here());
  }
}

}  // namespace ufo

