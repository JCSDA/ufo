/*
 * (C) Copyright 2017-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALTLAD_H_
#define UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALTLAD_H_

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variable.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/columnretrieval/ObsColumnRetrievalParameters.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// ColumnRetrieval TL/AD observation operator class
class ObsColumnRetrievalTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsColumnRetrievalTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsColumnRetrievalParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsColumnRetrievalTLAD";}

  ObsColumnRetrievalTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsColumnRetrievalTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;

  /// Copy the stored trajectory at location \p iobs into the per-location work arrays.
  void getTrajectoryAtLocation(size_t iobs, std::vector<double> & prsiModel,
                               std::vector<double> & avgKernel,
                               std::vector<double> & prsiObs) const;

  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  oops::Variable tracerVar_;
  oops::Variable prsiVar_;
  int nlayers_;
  bool isApriori_;
  bool isAveragingKernel_;
  bool totalNoVertice_;
  StretchVertices stretch_;
  double modelUnitsCoeff_;
  double tropopause_;

  // Trajectory storage (populated by setTrajectory)
  size_t nlocs_;
  size_t nval_;
  std::vector<double> avgKernelObs_;  // size nlayers_ * nlocs_
  std::vector<double> prsiObs_;       // size (nlayers_+1) * nlocs_
  std::vector<double> prsiModel_;     // size nval_ * nlocs_
  std::vector<bool> validLoc_;        // false where obs data are missing, size nlocs_
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALTLAD_H_
