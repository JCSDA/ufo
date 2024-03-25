/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_AEROSOLS_AODLUTS_OBSAODLUTS_H_
#define UFO_OPERATORS_AEROSOLS_AODLUTS_OBSAODLUTS_H_

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/Fortran.h"
#include "ufo/ObsOperatorBase.h"

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"

#include "ufo/operators/aerosols/AODLuts/ObsAodLUTs.interface.h"


namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

class AodLUTsOptionsParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(AodLUTsOptionsParameters, ObsOperatorParametersBase)
 public:
  oops::RequiredParameter<std::string> option{"AerosolOption", this};
  oops::RequiredParameter<std::string> file{"RCFile", this};
  oops::RequiredParameter<std::string> sensorID{"Sensor_ID", this};
  oops::OptionalParameter<std::string> coeffPath{"CoefficientPath", this};
  oops::OptionalParameter<std::string> endian{"EndianType", this};
  oops::Parameter<bool> absorptionAOD{"AbsorptionAod", false, this};
  oops::Parameter<double> modelunitscoeff{"model units coeff", 1, this};
  oops::Parameter<bool> drymixratio{"dry mix ratio", false, this};
};

class AodLUTsParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(AodLUTsParameters, ObsOperatorParametersBase)
 public:
  oops::RequiredParameter<AodLUTsOptionsParameters> options{"obs options", this};
};

// -----------------------------------------------------------------------------
/// AodLUTs observation for UFO.
class ObsAodLUTs : public ObsOperatorBase,
                    private util::ObjectCounter<ObsAodLUTs> {
 public:
  typedef AodLUTsParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsAodLUTs";}

  ObsAodLUTs(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodLUTs();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;


// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperAodLUTs_;}
  const int & toFortran() const {return keyOperAodLUTs_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperAodLUTs_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_AEROSOLS_AODLUTS_OBSAODLUTS_H_
