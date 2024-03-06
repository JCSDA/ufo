/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVAL_H_
#define UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVAL_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/columnretrieval/ObsColumnRetrieval.interface.h"
#include "ufo/operators/columnretrieval/ObsColumnRetrievalParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// ColumnRetrieval observation operator class
class ObsColumnRetrieval : public ObsOperatorBase,
                   private util::ObjectCounter<ObsColumnRetrieval> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsColumnRetrievalParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsColumnRetrieval";}

  ObsColumnRetrieval(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsColumnRetrieval();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVAL_H_
