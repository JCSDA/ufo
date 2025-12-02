/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>

#include "ufo/utils/dataextractor/DataExtractorInput.h"
#include "ufo/utils/dataextractor/DataExtractorNetCDFBackend.h"
#include "ufo/variabletransforms/Cal_SatRadianceFromPCScores.h"

#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace ufo {

/************************************************************************************/
//  Cal_SatRadianceFromPCScores
/************************************************************************************/

static TransformMaker<Cal_SatRadianceFromPCScores>
    makerCal_SatRadianceFromPCScores_("SatRadianceFromPCScores");

Cal_SatRadianceFromPCScores::Cal_SatRadianceFromPCScores(
    const Parameters_ &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr), parameters_(options),
      variables_({parameters_.pcVariable.value()}),
      channels_(parameters_.destinationVariable.value().channels()) {
  const size_t nComponents = parameters_.nComponents.value();
  for (std::size_t i = 0; i < nComponents; ++i) {
    const std::string str = std::to_string(i+1);
    variables_ += Variable(parameters_.pcVariable.value().fullName() + str);
  }
}

/************************************************************************************/

void Cal_SatRadianceFromPCScores::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> convert PC scores to spectral radiance"
            << std::endl;
  oops::Log::trace() << "      --> method: " << options_.Method.value() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  // Set up variables
  const size_t nlocs = obsdb_.nlocs();
  const size_t nComponents = parameters_.nComponents.value();

  // If no locations to process, return
  if (nlocs == 0) {
    return;
  }

  // Read reconstruction operator from netCDF file
  // NB the backend loadData method returns a 3-d array but operator dimensions are 2-d
  const DataExtractorNetCDFBackend<float> Extractor(parameters_.reconstructionFile.value());
  const std::string pcScoresGroup = parameters_.operatorDataGroup.value();
  const boost::multi_array<float, 3> payloadReconstructor
                                     = Extractor.loadData(pcScoresGroup).payloadArray;
  ASSERT(payloadReconstructor.shape()[2] == 1);  // 2-d data
  ASSERT(payloadReconstructor.shape()[0] == nComponents);  // no. PCs should match yaml
  const auto nReconstructorChannels = payloadReconstructor.shape()[1];
  Eigen::MatrixXf reconstructor
      = rowMajorMatrix::Map(payloadReconstructor.data(), nComponents, nReconstructorChannels);

  // Read channels associated with reconstruction operator from MetaData/sensorChannelNumber
  const std::string channelGroup = parameters_.operatorChannelGroup.value();
  const boost::multi_array<float, 3> payloadChannels
                                     = Extractor.loadData(channelGroup).payloadArray;
  ASSERT((payloadChannels.shape()[1] == 1) && (payloadChannels.shape()[2] == 1));  // 1-d data
  ASSERT(payloadChannels.num_elements() == nReconstructorChannels);
  std::vector<int> reconstructorChannels(nReconstructorChannels);
  for (std::size_t i = 0; i < nReconstructorChannels; ++i) {
    reconstructorChannels[i] = payloadChannels[i][0][0];
  }

  // If user specified Read Means associated with reconstruction operator from
  // MetaData/sensorChannelNumber
  std::vector<float> operatorMean(nReconstructorChannels);
  std::fill(operatorMean.begin(), operatorMean.end(), 0.0);

  if (parameters_.operatorMean.value().size() > 0) {
    const std::string meanGroup = parameters_.operatorMean.value();
    const boost::multi_array<float, 3> payloadMean
                                     = Extractor.loadData(meanGroup).payloadArray;
    ASSERT((payloadMean.shape()[1] == 1) && (payloadMean.shape()[2] == 1));  // 1-d data
    ASSERT(payloadMean.num_elements() == nReconstructorChannels);
    for (std::size_t i = 0; i < nReconstructorChannels; ++i) {
      operatorMean[i] = payloadMean[i][0][0];
    }
  }



  // Channels for destination will normally be the same as obsdb_.assimvariables().channels()
  const Variable radianceVar = parameters_.destinationVariable;
  const std::vector<int> destinationChannels = radianceVar.channels();
  const size_t nDestinationChannels = destinationChannels.size();

  // Find indices of destination channels within reconstruction operator channels
  std::vector<int> destinationChannelIndex(nReconstructorChannels);
  for (std::size_t i = 0; i < nDestinationChannels; ++i) {
    auto it = std::find(reconstructorChannels.begin(), reconstructorChannels.end(),
                        destinationChannels[i]);
    ASSERT(it != reconstructorChannels.end());
    destinationChannelIndex[i] = std::distance(reconstructorChannels.begin(), it);
  }

  // Read in principal component scores
  std::vector<std::vector<float>> pcScores(nComponents,
                                           std::vector<float>(nlocs, missingValueFloat));
  rowMajorMatrix pcScoresMatrix(nComponents, nlocs);
  for (std::size_t i = 0; i < nComponents; ++i) {
    std::string str = std::to_string(i+1);
    Variable pcvar = Variable(parameters_.pcVariable.value().fullName() + str);
    obsdb_.get_db(pcvar.group(), pcvar.variable(), pcScores[i]);
    pcScoresMatrix.row(i) = Eigen::Map<Eigen::VectorXf>(pcScores[i].data(), nlocs);
  }

  // Radiances computed as matrix multiplication of reconstruction operator and PC scores,
  // multiplying also by a scaling factor to produce radiances in W m^-2 sr^-1 (m^-1)^-1
  float scalingFactor = parameters_.scalingFactor.value();
  Eigen::MatrixXf reconMatrix \
    = reconstructor(Eigen::all, destinationChannelIndex).transpose() * pcScoresMatrix;
  reconMatrix *= scalingFactor;

  // Loop over channels and obs locations to populate derived observation variable
  std::vector<std::vector<float>> reconRadiance(nDestinationChannels,
                                                std::vector<float>(nlocs, missingValueFloat));
  for (size_t ichan = 0; ichan < nDestinationChannels; ++ichan) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      // Only for locations selected by where clause
      if (apply[iloc]) {
        reconRadiance[ichan][iloc] = reconMatrix(ichan, iloc) +\
                                     operatorMean[destinationChannelIndex[ichan]];
      }
    }
  }

  //  Write out the radiances to the Derived group and update qcflags
  for (size_t ichan = 0; ichan < nDestinationChannels; ++ichan) {
    putObservation(radianceVar.variable() + "_" + std::to_string(destinationChannels[ichan]),
                   reconRadiance[ichan]);
  }

  oops::Log::trace() << "Cal_SatRadianceFromPCScores::runTransform done" << std::endl;
}  // runTransform
}  // namespace ufo
