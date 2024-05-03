/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_SATTCWV_SATTCWVTLAD_H_
#define UFO_OPERATORS_SATTCWV_SATTCWVTLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// Precipitable water observation operator
class SatTCWVTLADParameters : public ObsOperatorParametersBase {
        OOPS_CONCRETE_PARAMETERS(SatTCWVTLADParameters, ObsOperatorParametersBase)
        // No additional option defined in YAML
};
class SatTCWVTLAD : public LinearObsOperatorBase,
                      private util::ObjectCounter<SatTCWVTLAD> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::SatTCWVTLAD";}

  typedef SatTCWVTLADParameters Parameters_;
  SatTCWVTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~SatTCWVTLAD();
  // Obs Operators
  // -----------------------------------------------------------------------------
  /*! \brief Compute Jacobian matrix d(TCWV)/d(q).
  *
  * \details This is based on Roger Saunders' Fortran original "sattcwv" 
  * observation operator code. This method must be called before calling the TL/AD
  * methods.
  *
  * \date Sept. 2021: Created by J. Hocking (Met Office)
  */
  // -----------------------------------------------------------------------------
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;

  // -----------------------------------------------------------------------------
  /*! \brief Given an increment to the model state, calculate an increment to the
  *   observation.
  *
  * \details This is based on Roger Saunders' Fortran original "sattcwv" 
  * observation operator code.
  *
  * \date Sept. 2021: Created by J. Hocking (Met Office)
  */
  // -----------------------------------------------------------------------------
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;

  // -----------------------------------------------------------------------------
  /*! \brief Given an increment to the observation, calculate the equivalent
  *   increment to the model state.
  *
  * \details This is based on Roger Saunders' Fortran original "sattcwv" 
  * observation operator code.
  *
  * \date Sept. 2021: Created by J. Hocking (Met Office)
  */
  // -----------------------------------------------------------------------------
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;



  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<const oops::Variables> varin_;

  std::vector<std::vector<double>> k_matrix;
  bool traj_init;
  size_t nlevels;
  size_t nprofiles;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_SATTCWV_SATTCWVTLAD_H_
