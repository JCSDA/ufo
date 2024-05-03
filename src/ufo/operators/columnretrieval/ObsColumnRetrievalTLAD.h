/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALTLAD_H_
#define UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALTLAD_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/columnretrieval/ObsColumnRetrievalParameters.h"
#include "ufo/operators/columnretrieval/ObsColumnRetrievalTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

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
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperColumnRetrieval_;}
  const int & toFortran() const {return keyOperColumnRetrieval_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperColumnRetrieval_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALTLAD_H_
