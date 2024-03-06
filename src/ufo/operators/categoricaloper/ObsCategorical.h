/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_CATEGORICALOPER_OBSCATEGORICAL_H_
#define UFO_OPERATORS_CATEGORICALOPER_OBSCATEGORICAL_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/categoricaloper/ObsCategoricalData.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Categorical observation operator.
///
/// The Categorical operator can be used to run several observation operators,
/// each of which produces a vector of H(x) values.
/// The Categorical operator then creates a final H(x) vector by selecting the
/// observation operator at each location according to a categorical variable.
///
/// The choice of observation operator at each location is governed by the `categorical variable`
/// parameter, which must be an integer or string variable in the MetaData group.
///
/// The `categorised operators` map is used to produce a correspondence between values of the
/// categorical variable and the operator used.
///
/// The `fallback operator` parameter governs the observation operator that will be used
/// whenever a particular value of the categorical variable does not exist in
/// `categorised operators`.
///
/// The `operator configurations` parameter governs the configuration of each of the operators
/// to be used.
/// If either the fallback operator or one of the categorised operators have not been configured,
/// an exception will be thrown.
///
/// An example yaml configuration is as follows:
///  obs operator:
///    name: Categorical
///    categorical variable: stationIdentification
///    fallback operator: "Composite"
///    categorised operators: {"47418": "Composite", "54857": "Identity"}
///    operator configurations:
///    - name: Identity
///    - name: Composite
///      components:
///       - name: Identity
///         variables:
///         - name: airTemperature
///         - name: surfacePressure
///       - name: VertInterp
///         variables:
///         - name: windNorthward
///         - name: windEastward
///
/// This operator uses MetaData/stationIdentification as the categorical variable.
/// Both the Identity and Composite operators are used to produce H(x) vectors.
/// Then, at each location in the ObsSpace:
/// - if MetaData/stationIdentification is equal to 47418 then the Composite H(x) is selected;
/// - if MetaData/stationIdentification is equal to 54857 then the Identity H(x) is selected;
/// - otherwise, the fallback operator (also Composite in this case) H(x) is selected.
class ObsCategorical : public ObsOperatorBase,
  private util::ObjectCounter<ObsCategorical> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsCategoricalParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsCategorical";}

  ObsCategorical(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsCategorical() override;

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return data_.requiredVars(); }

 private:
  void print(std::ostream &) const override;

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Data handler for the Categorical operator and TL/AD code.
  ObsCategoricalData<ObsOperatorBase> data_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_CATEGORICALOPER_OBSCATEGORICAL_H_
