/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_SATRADIANCEFROMPCSCORES_H_
#define UFO_VARIABLETRANSFORMS_CAL_SATRADIANCEFROMPCSCORES_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/Variable.h"
#include "ufo/filters/VariableTransformParametersBase.h"
#include "ufo/variabletransforms/TransformBase.h"



namespace ufo {

/// Configuration parameters for PC scores to radiance variable transformation
class Cal_SatRadianceFromPCScoresParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_SatRadianceFromPCScoresParameters,
                           VariableTransformParametersBase);

 public:
  /// The variable that will be created by the transform in the DerivedObsValue group.
  /// A list of channels to be simulated is specified. Example:
  ///
  /// destination variable:
  ///   name: radiance
  ///   channels: 1-20
  oops::RequiredParameter<Variable> destinationVariable{"destination variable", this};

  /// Path to netCDF file containing the principal component score reconstruction operator.
  /// It is expected that the operator is a two-dimensional variable in a group named (by default)
  /// PCScores with PC component as its first dimension and channel as its second dimension.
  /// There should be a single variable in a separate group named (by default) MetaData defining
  /// sensor channel number.
  /// Example contents:
  ///
  /// netcdf {
  /// dimensions:
  ///     Channel = 8461 ;
  ///     Component = 300 ;
  /// variables:
  ///     int Channel(Channel) ;
  ///     int Component(Component) ;
  /// group: MetaData {
  /// variables:
  ///     int sensorChannelNumber(Channel) ;
  /// } // group MetaData
  /// group: PCScores {
  /// variables:
  ///     float reconstructionOperator(Component, Channel) ;
  /// } // group PCScores
  /// }
  oops::RequiredParameter<std::string> reconstructionFile{"reconstruction operator file", this};

  /// Group in reconstruction operator file containing a 2-d variable which defines the operator,
  /// dimensioned by PC component and channel, by default this group is "PCScores"
  oops::Parameter<std::string> operatorDataGroup{"operator data group", "PCScores", this};

  /// Group in reconstruction operator file containing a 1-d variable defining channels,
  /// by default this group is "MetaData" and by convention the variable is "sensorChannelNumber"
  oops::Parameter<std::string> operatorChannelGroup{"operator channel group", "MetaData", this};

  /// Number of principal components to be used for reconstruction. This value should match
  /// both the number of components in the reconstruction operator file and the number of
  /// principal component scores variables in the obs space.
  oops::RequiredParameter<size_t> nComponents{"number of principal components", this};

  /// Variable containing PC scores. Example:
  ///
  /// principal component scores variable: MetaData/principalComponentScore
  ///
  /// In this case the obs space will contain MetaData/principalComponentScore1,
  /// MetaData/principalComponentScore2, ..., up to the number of principal components.
  oops::RequiredParameter<Variable> pcVariable{"principal component scores variable", this};

  /// Scaling factor applied to the matrix multiplication (reconstruction operator * PC scores)
  /// in order to reconstruct radiances in units of W m^-2 sr^-1 (m^-1)^-1.
  oops::Parameter<float> scalingFactor{"reconstruction scaling factor", 0.5, this};

  /// Group in reconstruction operator file containing a 1-d variable defining channels,
  /// by default this group is empty, unless the user specifies a PCMean group indicating
  /// that a mean radiance /bias should be added as part of the reconstruction.
  oops::Parameter<std::string> operatorMean{"operator mean data group", "", this};
};

/*!
* \brief Convert principal component scores to radiance.
*
* Reconstructs radiances by matrix multiplication of principal component scores
* with reconstruction operator.
*
* Example using all the parameter options:
*
* \code{.yaml}
* obs filters:
* - filter: Variable Transforms
*   Transform: SatRadianceFromPCScores
*   destination variable:
*     name: radiance
*     channels: *all_channels
*   reconstruction operator file: iasi_pcscores_reconstructor.nc4
*   operator data group: PCScores
*   operator channel group: MetaData
*   number of principal components: 300
*   principal component scores variable: MetaData/principalComponentScore
*   reconstruction scaling factor: 0.5
* \endcode
*
* See Cal_SatRadianceFromScaledRadianceParameters for filter setup.
*/

class Cal_SatRadianceFromPCScores : public TransformBase {
 public:
  typedef Cal_SatRadianceFromPCScoresParameters Parameters_;
  typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMajorMatrix;

  Cal_SatRadianceFromPCScores(const Parameters_ &options,
                              const ObsFilterData &data,
                              ioda::ObsDataVector<int> &flags,
                              ioda::ObsDataVector<float> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;
  Variables requiredVariables() const override { return variables_; }
 private:
  Parameters_ parameters_;
  Variables variables_;
  std::vector<int> channels_;
};

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_SATRADIANCEFROMPCSCORES_H_
