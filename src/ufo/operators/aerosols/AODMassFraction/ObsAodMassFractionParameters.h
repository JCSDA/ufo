/*
 *
 * Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_AEROSOLS_AODMASSFRACTION_OBSAODMASSFRACTIONPARAMETERS_H_
#define UFO_OPERATORS_AEROSOLS_AODMASSFRACTION_OBSAODMASSFRACTIONPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

enum class aerosolModel {
      UKCA_Dust
};

struct aerosolModelParameterTraitsHelper {
  typedef aerosolModel EnumType;
  static constexpr char enumTypeName[] = "aerosolModel";
  static constexpr util::NamedEnumerator<aerosolModel> namedValues[] = {
    { aerosolModel::UKCA_Dust, "UKCA_Dust" }
  };
};
}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::aerosolModel> :
    public EnumParameterTraits<ufo::aerosolModelParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

enum class extinctionMethod {
      MetOfficeLUTFit
};

struct extinctionMethodParameterTraitsHelper {
  typedef extinctionMethod EnumType;
  static constexpr char enumTypeName[] = "extinctionMethod";
  static constexpr util::NamedEnumerator<extinctionMethod> namedValues[] = {
    { extinctionMethod::MetOfficeLUTFit, "MetOfficeLUTFit" }
  };
};
}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::extinctionMethod> :
    public EnumParameterTraits<ufo::extinctionMethodParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// Configuration options recognized by the AodMassFraction operator.
class ObsAodMassFractionParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsAodMassFractionParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter<aerosolModel> aerosolModelName{
    "aerosolModel", "Aerosol model to be used", this};
  oops::RequiredParameter<extinctionMethod> extinctionMethodName{
    "extinctionMethod", "Method to calculate extinction coefficient", this};
  oops::OptionalParameter<std::vector<double>> coarseModeParams{
    "coarseModeParams", this};
  oops::OptionalParameter<std::vector<double>> accumulationModeParams{
    "accumulationModeParams", this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_AEROSOLS_AODMASSFRACTION_OBSAODMASSFRACTIONPARAMETERS_H_
