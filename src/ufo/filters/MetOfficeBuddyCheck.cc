/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/MetOfficeBuddyCheck.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <boost/format.hpp>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/sqr.h"
#include "ufo/filters/FilterUtils.h"
#include "ufo/filters/MetOfficeBuddyCheckParameters.h"
#include "ufo/filters/MetOfficeBuddyPair.h"
#include "ufo/filters/MetOfficeBuddyPairFinder.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"


namespace ufo {

namespace {

const double expArgMax = 80.0;  // maximum argument of the exponential function
const float maxGrossErrorProbability = 1.0;  // maximum allowed value for PGE

/// brief Extract 1st level data from the provided flattened vector.
template<typename T>
std::vector<T> extract1stLev(const std::vector<T> &flattenedArray,
                             const boost::optional<Eigen::ArrayXXi> &profileIndex) {
  if ( !profileIndex )
    return flattenedArray;

  std::vector<T> res;
  res.resize((*profileIndex).rows());
  for (size_t row=0; row < (*profileIndex).rows(); row++) {
    res[row] = flattenedArray[(*profileIndex)(row, 0)];
  }
  return res;
}

/// Return a 2D eigen array from the provided flattened vector.
Eigen::ArrayXXf unravel(const std::vector<float> &flattenedArray,
                        const boost::optional<Eigen::ArrayXXi> &profileIndex) {
  if (profileIndex) {
    Eigen::ArrayXXf res((*profileIndex).rows(), (*profileIndex).cols());
    for (size_t row=0; row < (*profileIndex).rows(); row++)
      for (size_t col=0; col < (*profileIndex).cols(); col++)
        res(row, col) = flattenedArray[(*profileIndex)(row, col)];
    return res;
  } else {
    Eigen::ArrayXXf res(flattenedArray.size(), 1);
    for (size_t row=0; row < flattenedArray.size(); row++)
      res(row, 0) = flattenedArray[row];
    return res;
  }
}

/// Update the provided flat array with values from the provided eigen array
void updateFlatData(std::vector<float> & flatArray, const Eigen::ArrayXXf &array,
                    const boost::optional<Eigen::ArrayXXi> &profileIndex) {
  if (profileIndex) {
    for (size_t row=0; row < (*profileIndex).rows(); row++)
      for (size_t col=0; col < (*profileIndex).cols(); col++)
        flatArray[(*profileIndex)(row, col)] = array(row, col);
  } else {
    ASSERT(array.cols() == 1);
    for (size_t row=0; row < array.rows(); row++)
      flatArray[row] = array(row, 0);
  }
}

/// \brief Given a vector of values associated with all unique locations held by all processes,
/// return a vector of values associated only with the locations held by the calling process.
///
/// \param globalValues
///   A vector of values associated with all unique locations held by all processes.
/// \param nlocs
///   Number of locations held by the calling process.
/// \param dist
///   An object managing the distribution of locations across processes.
std::vector<float> localValues(const std::vector<float> & globalValues,
                               size_t nlocs, const ioda::Distribution &dist) {
  std::vector<float> result(nlocs);
  for (size_t localObsId = 0; localObsId < nlocs; ++localObsId) {
    const size_t globalObsId = dist.globalUniqueConsecutiveLocationIndex(localObsId);
    result[localObsId] = globalValues[globalObsId];
  }
  return result;
}

/// Return a map mapping indices of records held on all MPI ranks to the vectors of global indices
/// of locations belonging to these records (sorted according to the criteria specified during
/// ObsSpace construction).
ioda::ObsSpace::RecIdxMap mapRecordIdsToLocations(const ioda::ObsSpace &obsdb,
                                                  const int numLevels,
                                                  const bool overrideObsGrouping) {
  const size_t nlocs = obsdb.nlocs();

  // Identify the index of each location within the record it belongs to
  std::vector<size_t> locationIndexWithinItsRecord(nlocs);
  for (ioda::ObsSpace::RecIdxIter irec = obsdb.recidx_begin(); irec != obsdb.recidx_end(); ++irec) {
    size_t index = 0;
    for (size_t loc : obsdb.recidx_vector(irec))
      locationIndexWithinItsRecord[loc] = index++;
  }

  // Identify the record containing each location
  std::vector<size_t> recordIds = obsdb.recnum();

  // Collect data from all MPI ranks
  obsdb.distribution()->allGatherv(locationIndexWithinItsRecord);
  obsdb.distribution()->allGatherv(recordIds);

  // Count the number of locations in each record
  std::map<size_t, size_t> numLocsPerRecord;
  for (size_t recordId : recordIds)
    ++numLocsPerRecord[recordId];

  // Map each record ID to the (ordered) indices of locations belonging to it
  ioda::ObsSpace::RecIdxMap result;

  if (overrideObsGrouping && numLevels == 1) {
    // Treat each location in each record separately in the buddy check even if
    // the ObsSpace has been divided into records.
    // This is only performed if the parameter `num_levels` has been set to 1.
    // The assignment of record ID works for all ObsSpace MPI distributions because at this
    // point in the code the `allGatherv` routine has been used to gather the obs
    // onto all processors.
    // In other words, recordIds[k] = k * obsdb.comm().size() for all PEs.
    for (size_t gloc = 0; gloc < recordIds.size(); ++gloc) {
      const size_t recid = gloc * obsdb.comm().size();
      result[recid] = {gloc};
    }
  } else {
    // Default behaviour
    for (const auto &kv : numLocsPerRecord) {
      const size_t recordId = kv.first;
      const size_t numLocsInRecord = kv.second;
      result[recordId].resize(numLocsInRecord);
    }
    for (size_t gloc = 0; gloc != recordIds.size(); ++gloc)
      result[recordIds[gloc]][locationIndexWithinItsRecord[gloc]] = gloc;
  }
  return result;
}

/// Return an array whose rows store the indices of locations belonging to successive profiles
/// that should be buddy-checked.
///
/// The buddy check operates on all profiles unless the ObsSpace contains the
/// `MetaData/extendedObsSpace` variable, in which case only profiles for which this variable is
/// set to 1 are buddy-checked. All profiles to be buddy-checked must comprise exactly \p numLevels
/// locations; an exception is thrown if that is not the case.
///
Eigen::ArrayXXi deriveIndices(const ioda::ObsSpace & obsdb,
                              const int numLevels,
                              const bool overrideObsGrouping) {
  // Assume ObsSpace contains only the averaged profiles if this variable isn't present.
  boost::optional<std::vector<int>> extended_obs_space;
  if (obsdb.has("MetaData", "extendedObsSpace")) {
    extended_obs_space = std::vector<int>(obsdb.nlocs());
    obsdb.get_db("MetaData", "extendedObsSpace", *extended_obs_space);
    obsdb.distribution()->allGatherv(*extended_obs_space);
  }

  ioda::ObsSpace::RecIdxMap locationsPerRecord =
    mapRecordIdsToLocations(obsdb, numLevels, overrideObsGrouping);
  Eigen::ArrayXXi profileIndex{locationsPerRecord.size(), numLevels};

  int recnum = 0;
  for (const auto &recordIndexAndLocations : locationsPerRecord) {
    const std::vector<size_t> &locations = recordIndexAndLocations.second;
    int levnum = 0;
    for (size_t loc : locations) {
      if (extended_obs_space && ((*extended_obs_space)[loc] != 1))
        continue;
      profileIndex(recnum, levnum) = loc;
      levnum++;
    }
    if (levnum == numLevels) {
      recnum++;
    } else if (levnum == 0) {
      // Ignore this profile
    } else {
      std::stringstream msg;
      msg << "Record (profile): " << recnum << " length: " << levnum << " does not match the "
          << "number of levels expected: " << numLevels;
      throw eckit::UserError(msg.str(), Here());
    }
  }
  profileIndex.conservativeResize(recnum, Eigen::NoChange);
  return profileIndex;
}

std::string fullVariableName(const Variable &var)
{
  if (var.group().empty())
    return var.variable();
  else
    return var.group() + "/" + var.variable();
}

template <typename T>
std::vector<T> uniqueElements(const std::vector<T> &v)
{
  std::vector<T> result(v);
  std::sort(result.begin(), result.end());
  const auto newEnd = std::unique(result.begin(), result.end());
  result.erase(newEnd, result.end());
  return result;
}

template <typename T>
std::vector<int> mapDistinctValuesToDistinctInts(const std::vector<T> &values)
{
  const std::vector<T> uniqueValues = uniqueElements(values);
  std::map<T, int> intByValue;
  {
    int index = 0;
    for (const T& value : uniqueValues)
      intByValue[value] = index++;
  }

  std::vector<int> ints;
  ints.reserve(values.size());
  std::transform(values.begin(), values.end(), std::back_inserter(ints),
                 [&](const T& value) { return intByValue.at(value); });
  return ints;
}

/// Data associated with a single filter variable.
struct ScalarVariableData {
  // If the buddy check operates on individual locations, the following arrays have dimensions
  // (`gnlocs`, 1), where `gnlocs` is the total number of unique locations held by all processes.
  // If the buddy check operates on profiles (e.g. radiosonde soundings), the arrays have
  // dimensions (`gnprofs`, `nlevels`), where `gnprofs` is the total number of buddy-checked
  // profiles held by all processes and `nlevels` is the number of levels per profile (taken from
  // the YAML option `num_levels`).
  Eigen::ArrayXXf obsValues;
  Eigen::ArrayXXf obsErrors;
  Eigen::ArrayXXf bgValues;
  Eigen::ArrayXXf bgErrors;
  Eigen::ArrayXXf bgErrors2;
  Eigen::ArrayXXf grossErrorProbabilities;

  // If the buddy check operates on individual locations, this vector has length `gnlocs` and
  // stores the QC flags at all locations. Otherwise it has length `gnprofs` and stores the
  // QC flags at the first location in each profile.
  std::vector<int> varFlags;

  // Gross error probabilities at all unique locations held by any process.
  std::vector<float> flatGrossErrorProbabilities;
};


}  // namespace

/// \brief Metadata of all observations processed by the filter.
struct MetOfficeBuddyCheck::MetaData {
  std::vector<float> latitudes;
  std::vector<float> longitudes;
  std::vector<util::DateTime> datetimes;
  boost::optional<std::vector<float>> pressures;
  boost::optional<Eigen::ArrayXXf> pressuresML;
  std::vector<int> stationIds;
};

MetOfficeBuddyCheck::MetOfficeBuddyCheck(ioda::ObsSpace& obsdb, const Parameters_& parameters,
                                         std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                         std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, std::move(flags), std::move(obserr),
                VariableNameMap(parameters.AliasFile.value())), options_(parameters)
{
  oops::Log::debug() << "MetOfficeBuddyCheck: config = " << options_ << std::endl;
  allvars_ += Variables(filtervars_, "HofX");
  for (size_t i = 0; i < filtervars_.size(); ++i) {
    allvars_ += backgroundErrorVariable(filtervars_[i],
                                        options_.backgroundErrorSuffix.value(),
                                        options_.backgroundErrorGroup.value());
    if (options_.backgroundErrorGroup2.value()) {
      allvars_ += backgroundErrorVariable(filtervars_[i],
                                          options_.backgroundErrorSuffix2.value(),
                                          options_.backgroundErrorGroup2.value().value());
    }
  }
}

void MetOfficeBuddyCheck::applyFilter(const std::vector<bool> & apply,
                                      const Variables & filtervars,
                                      std::vector<std::vector<bool>> & flagged) const {
  // Fetch metadata required for identifying buddy pairs.
  const boost::optional<int> &numLevels = options_.numLevels.value();
  const bool overrideObsGrouping = options_.overrideObsGrouping;
  boost::optional<Eigen::ArrayXXi> profileIndex;
  if (numLevels)
    profileIndex = deriveIndices(obsdb_, *numLevels, overrideObsGrouping);

  const std::vector<size_t> validObsIds = getValidObservationIds(apply, profileIndex);
  MetaData obsData = collectMetaData(profileIndex);

  const std::vector<float> bgErrorHorizCorrScales = calculatePiecewiseLinear(
                                validObsIds,
                                obsData.latitudes,
                                options_.horizontalCorrelationScaleInterpolationPoints.value());
  std::vector<float> bgErrorHorizCorrScales2;
  std::vector<float> anisotropy;
  std::vector<float> anisotropy2;
  if (options_.horizontalCorrelationScale2InterpolationPoints.value()) {
    bgErrorHorizCorrScales2 = calculatePiecewiseLinear(
                          validObsIds,
                          obsData.latitudes,
                          options_.horizontalCorrelationScale2InterpolationPoints.value().value());
  }
  if (options_.anisotropyInterpolationPoints.value()) {
    anisotropy = calculatePiecewiseLinear(
                          validObsIds,
                          obsData.latitudes,
                          options_.anisotropyInterpolationPoints.value().value());
  }
  if (options_.anisotropy2InterpolationPoints.value()) {
    anisotropy2 = calculatePiecewiseLinear(
                          validObsIds,
                          obsData.latitudes,
                          options_.anisotropy2InterpolationPoints.value().value());
  }

  const std::vector<bool> verbose = flagAndPrintVerboseObservations(
        validObsIds, obsData.latitudes, obsData.longitudes, obsData.datetimes,
        obsData.pressures.get_ptr(), obsData.stationIds, bgErrorHorizCorrScales);

  // Identify buddy pairs

  const std::vector<float> *pressures =
      options_.sortByPressure ? obsData.pressures.get_ptr() : nullptr;
  MetOfficeBuddyPairFinder buddyPairFinder(options_, obsData.latitudes, obsData.longitudes,
                                           obsData.datetimes, pressures,
                                           obsData.stationIds);
  const std::vector<MetOfficeBuddyPair> buddyPairs = buddyPairFinder.findBuddyPairs(validObsIds);

  std::shared_ptr<const ioda::Distribution> distribution = obsdb_.distribution();

  // Fetch data/metadata required for buddy-check calculation.

  auto getFilterVariableName = [&] (size_t filterVarIndex) {
    return filtervars.variable(filterVarIndex).variable();
  };

  // Return a ScalarVariableData object containing the data (observed and background values, errors
  // etc.) corresponding to a specific filter variable.
  auto getScalarVariableData = [&] (size_t filterVarIndex) {
    // Get the name of the filter variable.
    const std::string filterVar = getFilterVariableName(filterVarIndex);

    // Retrieve data.
    const std::vector<float> obsValues =
        getGlobalVariable<float>("ObsValue", filterVar);
    const std::vector<float> obsErrors =
        getGlobalVariable<float>("ObsErrorData", filterVar);
    const std::vector<float> bgValues =
        getGlobalVariable<float>("HofX", filterVar);
    const std::vector<float> bgErrors =
        getGlobalVariable<float>(backgroundErrorVariable(filtervars[filterVarIndex],
                                                         options_.backgroundErrorSuffix.value(),
                                                         options_.backgroundErrorGroup.value()));
    std::vector<float> bgErrors2(bgErrors.size(), util::missingValue<float>());
    if (options_.backgroundErrorGroup2.value()) {
      bgErrors2 = getGlobalVariable<float>(backgroundErrorVariable(filtervars[filterVarIndex],
                                              options_.backgroundErrorSuffix2.value(),
                                              options_.backgroundErrorGroup2.value().value()));
    }
    const std::vector<int> flags =
        getGlobalVariable<int>("QCflagsData", filterVar);
    const std::vector<float> grossErrorProbabilities =
        getGlobalVariable<float>("GrossErrorProbability", filterVar);

    // Store data in a ScalarVariableData object.
    ScalarVariableData data;

    // If we're buddy-checking profiles, unravel vectors into 2D arrays indexed by profile.
    data.obsValues = unravel(obsValues, profileIndex);
    data.obsErrors = unravel(obsErrors, profileIndex);
    data.bgValues = unravel(bgValues, profileIndex);
    data.bgErrors = unravel(bgErrors, profileIndex);
    data.bgErrors2 = unravel(bgErrors2, profileIndex);
    data.grossErrorProbabilities = unravel(grossErrorProbabilities, profileIndex);

    data.varFlags = extract1stLev(flags, profileIndex);

    data.flatGrossErrorProbabilities = std::move(grossErrorProbabilities);

    return data;
  };

  // Buddy-check all filter variables

  // Gross error probabilities updated by buddy check, indexed by variable name.
  // (The updated values are copied to the ObsSpace only after all variables have been processed).
  std::map<std::string, std::vector<float>> calculatedGrossErrProbsByVarName;

  bool previousVariableWasFirstComponentOfTwo = false;
  for (size_t filterVarIndex = 0; filterVarIndex < filtervars.size(); ++filterVarIndex) {
    if (previousVariableWasFirstComponentOfTwo) {
      // Vector (two-component) variable
      ScalarVariableData firstComponentData =
          getScalarVariableData(filterVarIndex - 1);
      ScalarVariableData secondComponentData =
          getScalarVariableData(filterVarIndex);

      for (Eigen::Index i = 0; i < firstComponentData.grossErrorProbabilities.rows(); ++i) {
        for (Eigen::Index j = 0; j < firstComponentData.grossErrorProbabilities.cols(); ++j) {
          firstComponentData.grossErrorProbabilities(i, j) =
              std::max(firstComponentData.grossErrorProbabilities(i, j),
                       secondComponentData.grossErrorProbabilities(i, j));
        }
      }

      checkVectorData(buddyPairs, firstComponentData.varFlags, verbose,
                      bgErrorHorizCorrScales, obsData.stationIds, obsData.datetimes,
                      obsData.pressuresML.get_ptr(), firstComponentData.obsValues,
                      secondComponentData.obsValues, firstComponentData.obsErrors,
                      firstComponentData.bgValues, secondComponentData.bgValues,
                      firstComponentData.bgErrors, firstComponentData.grossErrorProbabilities);
      // OPS doesn't update the gross error probabilities of the second component variable,
      // but it seems more consistent to do so (and it facilitates the implementation of
      // flagRejectedObservations().
      // The implementation of checkVectorSurfaceData() still assumes that the *input* gross error
      // probabilities, flags and background error estimates are the same for both components.

      // Update the flat PGEs and extract those associated with locations held by this process.
      updateFlatData(firstComponentData.flatGrossErrorProbabilities,
                     firstComponentData.grossErrorProbabilities, profileIndex);
      std::vector<float> localFlatGrossErrorProbabilities = localValues(
            firstComponentData.flatGrossErrorProbabilities, obsdb_.nlocs(), *distribution);
      calculatedGrossErrProbsByVarName[getFilterVariableName(filterVarIndex - 1)] =
          localFlatGrossErrorProbabilities;
      calculatedGrossErrProbsByVarName[getFilterVariableName(filterVarIndex)] =
          std::move(localFlatGrossErrorProbabilities);

      previousVariableWasFirstComponentOfTwo = false;
    } else {
      if (filtervars[filterVarIndex].options().getBool("first_component_of_two", false)) {
        previousVariableWasFirstComponentOfTwo = true;
      } else {
        // Scalar variable
        ScalarVariableData data =
            getScalarVariableData(filterVarIndex);
        checkScalarData(buddyPairs, data.varFlags, verbose, bgErrorHorizCorrScales,
                        options_.horizontalCorrelationScale2InterpolationPoints.value()?
                        &bgErrorHorizCorrScales2 : nullptr,
                        options_.anisotropyInterpolationPoints.value()?
                        &anisotropy : nullptr,
                        options_.anisotropy2InterpolationPoints.value()?
                        &anisotropy2 : nullptr,
                        obsData.stationIds, obsData.datetimes, obsData.pressuresML.get_ptr(),
                        data.obsValues, data.obsErrors,
                        data.bgValues, data.bgErrors,
                        options_.backgroundErrorGroup2.value()?
                        &data.bgErrors2 : nullptr,
                        data.grossErrorProbabilities);

        // Update the flat PGEs and extract those associated with locations held by this process.
        updateFlatData(data.flatGrossErrorProbabilities,
                       data.grossErrorProbabilities, profileIndex);
        std::vector<float> localFlatGrossErrorProbabilities = localValues(
              data.flatGrossErrorProbabilities, obsdb_.nlocs(), *distribution);
        calculatedGrossErrProbsByVarName[getFilterVariableName(filterVarIndex)] =
            std::move(localFlatGrossErrorProbabilities);
      }
    }
  }

  // Update observations flags and gross error probabilities

  for (const auto &varNameAndGrossErrProbs : calculatedGrossErrProbsByVarName)
    obsdb_.put_db("GrossErrorProbability", varNameAndGrossErrProbs.first,
                  varNameAndGrossErrProbs.second);

  flagRejectedObservations(filtervars, calculatedGrossErrProbsByVarName, flagged);
}

Variable MetOfficeBuddyCheck::backgroundErrorVariable(const Variable &filterVariable,
                                                      const std::string &suffix,
                                                      const std::string &groupName) const {
  oops::Log::debug() << "BGE var: " << groupName + "/" +
                        nameMap_.convertName(filterVariable.variable()).name() +
                        suffix << std::endl;
  return Variable(groupName + "/" + nameMap_.convertName(filterVariable.variable()).name()
                  + suffix);
}

MetOfficeBuddyCheck::MetaData MetOfficeBuddyCheck::collectMetaData(
    const boost::optional<Eigen::ArrayXXi> & profileIndex) const {
  MetaData obsData;

  std::shared_ptr<const ioda::Distribution> distribution = obsdb_.distribution();

  obsData.latitudes = getGlobalObsSpaceVariable<float>("MetaData", "latitude");
  obsData.longitudes = getGlobalObsSpaceVariable<float>("MetaData", "longitude");
  obsData.datetimes = getGlobalObsSpaceVariable<util::DateTime>("MetaData", "dateTime");

  if (obsdb_.has(options_.pressureGroup, options_.pressureCoord)) {
    obsData.pressures =
      getGlobalObsSpaceVariable<float>(options_.pressureGroup, options_.pressureCoord);
    obsData.pressuresML = unravel(*obsData.pressures, profileIndex);
  }
  obsData.stationIds = getStationIds();

  if (profileIndex) {
    obsData.latitudes = extract1stLev(obsData.latitudes, profileIndex);
    obsData.longitudes = extract1stLev(obsData.longitudes, profileIndex);
    obsData.datetimes = extract1stLev(obsData.datetimes, profileIndex);
    obsData.stationIds = extract1stLev(obsData.stationIds, profileIndex);
    if (obsdb_.has(options_.pressureGroup, options_.pressureCoord)) {
      obsData.pressures = extract1stLev(*obsData.pressures, profileIndex);
    }
  }
  return obsData;
}

std::vector<int> MetOfficeBuddyCheck::getStationIds() const {
  const boost::optional<Variable> &stationIdVariable = options_.stationIdVariable.value();
  if (stationIdVariable == boost::none) {
    std::vector<int> stationIds;
    if (obsdb_.obs_group_vars().empty()) {
      // Observations were not grouped into records.
      // Assume all observations were taken by the same station.
      stationIds.assign(obsdb_.nlocs(), 0);
    } else {
      const std::vector<size_t> &recordNumbers = obsdb_.recnum();
      stationIds.assign(recordNumbers.begin(), recordNumbers.end());
    }
    obsdb_.distribution()->allGatherv(stationIds);
    return stationIds;
  } else {
    switch (obsdb_.dtype(stationIdVariable->group(), stationIdVariable->variable())) {
    case ioda::ObsDtype::Integer:
      return getGlobalVariable<int>(stationIdVariable->group(), stationIdVariable->variable());

    case ioda::ObsDtype::String:
      {
        std::vector<std::string> stringIds = getGlobalVariable<std::string>(
              stationIdVariable->group(), stationIdVariable->variable());
        return mapDistinctValuesToDistinctInts(stringIds);
      }

    default:
      throw eckit::UserError("Only integer and string variables may be used as station IDs",
                             Here());
    }
  }
}

template <typename T>
std::vector<T> MetOfficeBuddyCheck::getGlobalObsSpaceVariable(const std::string &group,
                                                              const std::string &variable) const {
  std::vector<T> values(obsdb_.nlocs());
  obsdb_.get_db(group, variable, values);
  obsdb_.distribution()->allGatherv(values);
  return values;
}

template <typename T>
std::vector<T> MetOfficeBuddyCheck::getGlobalVariable(const std::string &group,
                                                      const std::string &variable) const {
  return getGlobalVariable<T>(Variable(group + "/" + variable));
}

template <typename T>
std::vector<T> MetOfficeBuddyCheck::getGlobalVariable(const Variable &var) const {
  std::vector<T> values;
  data_.get(var, values);
  obsdb_.distribution()->allGatherv(values);
  return values;
}

std::vector<float> MetOfficeBuddyCheck::calculatePiecewiseLinear(
    const std::vector<size_t> &validObsIds,
    const std::vector<float> &latitudes,
    const std::map<float, float> &interpolationPoints) const {

  std::vector<double> abscissas, ordinates;
  for (const std::pair<const float, float> &xy :
       interpolationPoints) {
    abscissas.push_back(xy.first);
    ordinates.push_back(xy.second);
  }
  PiecewiseLinearInterpolation horizontalCorrelationScaleInterp(std::move(abscissas),
                                                                std::move(ordinates));

  std::vector<float> scales(latitudes.size());
  for (size_t obsId : validObsIds) {
      scales[obsId] = horizontalCorrelationScaleInterp(latitudes[obsId]);
  }

  return scales;
}

std::vector<bool> MetOfficeBuddyCheck::flagAndPrintVerboseObservations(
    const std::vector<size_t> &validObsIds,
    const std::vector<float> &latitudes,
    const std::vector<float> &longitudes,
    const std::vector<util::DateTime> &datetimes,
    const std::vector<float> *pressures,
    const std::vector<int> &stationIds,
    const std::vector<float> &bgErrorHorizCorrScales) const {

  size_t numVerboseObs = 0;
  std::vector<bool> verbose(latitudes.size(), false);

  for (size_t obsId : validObsIds) {
    verbose[obsId] = std::any_of(
          options_.tracedBoxes.value().begin(),
          options_.tracedBoxes.value().end(),
          [&](const LatLonBoxParameters &box)
          { return box.contains(latitudes[obsId], longitudes[obsId]); });

    if (verbose[obsId]) {
      if (numVerboseObs == 0) {
        oops::Log::trace() << "Obs   StationId Lat     Lon     Pressure "
                              "Datetime             CorrScale\n";
      }

      ++numVerboseObs;
      const float pressure = pressures != nullptr ? (*pressures)[obsId] : 0;
      oops::Log::trace() << boost::format("%5d %9d %7.2f %7.2f %8.0f %s %9.1f\n") %
                            obsId % stationIds[obsId] % latitudes[obsId] % longitudes[obsId] %
                            pressure % datetimes[obsId] % bgErrorHorizCorrScales[obsId];
    }
  }

  oops::Log::trace() << "Buddy check: " << numVerboseObs << " verbose observations" << std::endl;

  return verbose;
}

void MetOfficeBuddyCheck::checkScalarData(const std::vector<MetOfficeBuddyPair> &pairs,
                                          const std::vector<int> &flags,
                                          const std::vector<bool> &verbose,
                                          const std::vector<float> &bgErrorHorizCorrScales,
                                          const std::vector<float> *bgErrorHorizCorrScales2,
                                          const std::vector<float> *anisotropy,
                                          const std::vector<float> *anisotropy2,
                                          const std::vector<int> &stationIds,
                                          const std::vector<util::DateTime> &datetimes,
                                          const Eigen::ArrayXXf *pressures,
                                          const Eigen::ArrayXXf &obsValues,
                                          const Eigen::ArrayXXf &obsErrors,
                                          const Eigen::ArrayXXf &bgValues,
                                          const Eigen::ArrayXXf &bgErrors,
                                          const Eigen::ArrayXXf *bgErrors2,
                                          Eigen::ArrayXXf &pges) const {
  using util::sqr;
  const boost::optional<int> &nolevs = options_.numLevels.value();
  const int numLevels = nolevs ? *nolevs : 0;  // number of actual levels
  const Eigen::Index cols = obsValues.cols();  // number of data columns (levels) to loop over

  const bool isMaster = obsdb_.comm().rank() == 0;
  if (isMaster) {
    oops::Log::trace() << "ObsA  ObsB  StatIdA  StatIdB  lev  DiffA DiffB "
                          "Dist   Corr  Agree   PgeA   PgeB   Mult\n";
  }

  const double invTemporalCorrScale = 1.0 / options_.temporalCorrelationScale.value().toSeconds();
  const float missing = util::missingValue<float>();

  // Loop over buddy pairs
  for (const MetOfficeBuddyPair &pair : pairs) {
    const size_t jA = pair.obsIdA;
    const size_t jB = pair.obsIdB;

    if (bgErrorHorizCorrScales2 && anisotropy && anisotropy2 && bgErrors2) {  // 2 length scales
      // - hcScale1: synoptic horizontal error scale for the pair of obs
      // - hcScale2: mesoscale horizontal error scale for the pair of obs
      // - sin2: sin-squared of angle between E-W and line joining pair of obs
      // - aniso1/2: synoptic/mesoscale anisotropy
      // - scaledDist1/2: synoptic/mesoscale scaled distance between observations
      const double hcScale1 = 0.5 * (bgErrorHorizCorrScales[jA] +
                                     bgErrorHorizCorrScales[jB]);
      const double hcScale2 = 0.5 * ((*bgErrorHorizCorrScales2)[jA] +
                                     (*bgErrorHorizCorrScales2)[jB]);
      const double sin2 = sqr(std::sin(0.5*(pair.rotationAInRad + pair.rotationBInRad)));
      const double aniso1 = 0.5 * ((*anisotropy)[jA] + (*anisotropy)[jB]);
      const double aniso2 = 0.5 * ((*anisotropy2)[jA] + (*anisotropy2)[jB]);
      const double scaledDist1 = std::pow(1 + sin2*(sqr(aniso1) - 1), 0.5) *
                                 pair.distanceInKm / hcScale1;
      const double scaledDist2 = std::pow(1 + sin2*(sqr(aniso2) - 1), 0.5) *
                                 pair.distanceInKm / hcScale2;

      // Background error correlation and covariance between observation positions.
      double covar;

      // Loop over each level
      for (Eigen::Index jlev=0; jlev < cols; jlev++) {
        if (pges(jA, jlev) >= maxGrossErrorProbability ||
            pges(jB, jlev) >= maxGrossErrorProbability ||
            pges(jA, jlev) == missing ||
            pges(jB, jlev) == missing ||
            obsValues(jA, jlev) == missing ||
            obsValues(jB, jlev) == missing ||
            bgValues(jA, jlev) == missing ||
            bgValues(jB, jlev) == missing ||
            bgErrors(jA, jlev) == missing ||
            bgErrors(jB, jlev) == missing ||
            (*bgErrors2)(jA, jlev) == missing ||
            (*bgErrors2)(jB, jlev) == missing
            )
          continue;  // skip to next level

        if (numLevels == 1) {
          // Single level data.
          covar = bgErrors(jA, jlev) * bgErrors(jB, jlev) *
                    (1.0 + scaledDist1) * std::exp(-scaledDist1 -
                    options_.verticalCorrelationScale.value() *
                    sqr(std::log((*pressures)(jA, 0)/(*pressures)(jB, 0))) -
                    sqr((datetimes[jA] - datetimes[jB]).toSeconds() * invTemporalCorrScale)) +
                  (*bgErrors2)(jA, jlev) * (*bgErrors2)(jB, jlev) *
                    (1.0 + scaledDist2) * std::exp(-scaledDist2 -
                    options_.verticalCorrelationScale.value() *
                    sqr(std::log((*pressures)(jA, 0)/(*pressures)(jB, 0))) -
                    sqr((datetimes[jA] - datetimes[jB]).toSeconds() * invTemporalCorrScale));
        } else {
          // Multi-level/surface data; treat vertical correlation as 1.0
          covar = bgErrors(jA, jlev) * bgErrors(jB, jlev) *
                    (1.0 + scaledDist1) * std::exp(-scaledDist1 -
                    sqr((datetimes[jA] - datetimes[jB]).toSeconds() * invTemporalCorrScale)) +
                  (*bgErrors2)(jA, jlev) * (*bgErrors2)(jB, jlev) *
                    (1.0 + scaledDist2) * std::exp(-scaledDist2 -
                    sqr((datetimes[jA] - datetimes[jB]).toSeconds() * invTemporalCorrScale));
        }
        // Differences from background
        const double diffA = obsValues(jA, jlev) - bgValues(jA, jlev);
        const double diffB = obsValues(jB, jlev) - bgValues(jB, jlev);
        // Estimated error variances (ob+bk) (eqn 2.5)
        const double errVarA = sqr(obsErrors(jA, jlev)) +
                               sqr(bgErrors(jA, jlev)) + sqr((*bgErrors2)(jA, jlev));
        const double errVarB = sqr(obsErrors(jB, jlev)) +
                               sqr(bgErrors(jB, jlev)) + sqr((*bgErrors2)(jB, jlev));
        // (Total error correlation between ob positions)**2
        const double rho2 = sqr(covar) / (errVarA * errVarB);
        // Argument for exponents
        double expArg = -(0.5 * rho2 / (1.0 - rho2)) *
            (sqr(diffA) / errVarA + sqr(diffB) / errVarB - 2.0 * diffA * diffB / covar);
        expArg = options_.dampingFactor1 * (-0.5 * std::log(1.0 - rho2) + expArg);  //
        expArg = std::min(expArgMax, std::max(-expArgMax, expArg));                 //
        // Z = P(OA)*P(OB)/P(OA and OB)
        double z = 1.0 / (1.0 - (1.0 - pges(jA, jlev)) * (1.0 - pges(jB, jlev)) *
            (1.0 - std::exp(expArg)));
        if (z <= 0.0)
          z = 1.0;  // rounding error control
        z = std::pow(z, options_.dampingFactor2);
        pges(jA, jlev) *= z;                       // update PGE
        pges(jB, jlev) *= z;                       // update PGE
        if (isMaster && (verbose[jA] || verbose[jB])) {
          const double corr = covar / std::pow(
                              (sqr(bgErrors(jA, jlev)) + sqr((*bgErrors2)(jA, jlev))) *
                              (sqr(bgErrors(jB, jlev)) + sqr((*bgErrors2)(jB, jlev))),
                              0.5);

          oops::Log::trace() << boost::format("%5d %5d %8d %8d %5d "
                                              "%5.3f %5.3f %6.3f "
                                              "%6.4f %6.3f %6.4f %6.4f %6.4f\n") %
                                jA % jB % stationIds[jA] % stationIds[jB] % jlev %
                                diffA % diffB % pair.distanceInKm %
                                corr % std::exp(expArg) % pges(jA, jlev) % pges(jB, jlev) % z;
        }  // verbose
      }  // jlev
    } else {  // single length scale
      if (!options_.useAllObservations) {
        // Check that observations are valid and buddy check is required
        if (flags[jA] != QCflags::pass || flags[jB] != QCflags::pass)
          continue;  // skip to next pair
      }
      // eqn 3.9
      // - hcScale: horizontal error scale for the pair of obs
      // - scaledDist: scaled Distance between observations
      const double hcScale = 0.5 * (bgErrorHorizCorrScales[jA] + bgErrorHorizCorrScales[jB]);
      const double scaledDist = pair.distanceInKm / hcScale;

      // Background error correlation between observation positions.
      // eqns 3.10, 3.11
      double corr;
      if (numLevels == 1) {
        // Single level data.
        corr = (1.0 + scaledDist) *
              std::exp(-scaledDist -
                        options_.verticalCorrelationScale.value() *
                        sqr(std::log((*pressures)(jA, 0)/(*pressures)(jB, 0))) -
                        sqr((datetimes[jA] - datetimes[jB]).toSeconds() * invTemporalCorrScale));
      } else {
        // Multi-level/surface data; treat vertical correlation as 1.0
        corr = (1.0 + scaledDist) *
              std::exp(-scaledDist -
                        sqr((datetimes[jA] - datetimes[jB]).toSeconds() * invTemporalCorrScale));
      }
      if (corr < 0.1)  // Check against minimum background error correlation.
        continue;  // skip to next pair

      // Loop over each level
      for (Eigen::Index jlev=0; jlev < cols; jlev++) {
        if (pges(jA, jlev) >= maxGrossErrorProbability ||
            pges(jB, jlev) >= maxGrossErrorProbability ||
            pges(jA, jlev) == missing ||
            pges(jB, jlev) == missing ||
            obsValues(jA, jlev) == missing ||
            obsValues(jB, jlev) == missing ||
            bgValues(jA, jlev) == missing ||
            bgValues(jB, jlev) == missing ||
            bgErrors(jA, jlev) == missing ||
            bgErrors(jB, jlev) == missing
            )
          continue;  // skip to next level

        // Differences from background
        const double diffA = obsValues(jA, jlev) - bgValues(jA, jlev);
        const double diffB = obsValues(jB, jlev) - bgValues(jB, jlev);
        // Estimated error variances (ob+bk) (eqn 2.5)
        const double errVarA = sqr(obsErrors(jA, jlev)) + sqr(bgErrors(jA, jlev));
        const double errVarB = sqr(obsErrors(jB, jlev)) + sqr(bgErrors(jB, jlev));
        // Background error covariance between ob positions (eqn 3.13)
        const double covar = corr * bgErrors(jA, jlev) * bgErrors(jB, jlev);
        // (Total error correlation between ob positions)**2 (eqn 3.14)
        const double rho2 = sqr(covar) / (errVarA * errVarB);
        // Argument for exponents
        double expArg = -(0.5 * rho2 / (1.0 - rho2)) *
            (sqr(diffA) / errVarA + sqr(diffB) / errVarB - 2.0 * diffA * diffB / covar);
        expArg = options_.dampingFactor1 * (-0.5 * std::log(1.0 - rho2) + expArg);  // exponent of
        expArg = std::min(expArgMax, std::max(-expArgMax, expArg));                 // eqn 3.18
        // Z = P(OA)*P(OB)/P(OA and OB)
        double z = 1.0 / (1.0 - (1.0 - pges(jA, jlev)) * (1.0 - pges(jB, jlev)) *
            (1.0 - std::exp(expArg)));
        if (z <= 0.0)
          z = 1.0;  // rounding error control
        z = std::pow(z, options_.dampingFactor2);  // eqn 3.16
        pges(jA, jlev) *= z;                       // eqn 3.17
        pges(jB, jlev) *= z;                       // eqn 3.17
        if (isMaster && (verbose[jA] || verbose[jB])) {
          oops::Log::trace() << boost::format("%5d %5d %8d %8d %5d "
                                              "%5.1f %5.1f %6.1f "
                                              "%5.3f %6.3f %6.3f %6.3f %6.3f\n") %
                                jA % jB % stationIds[jA] % stationIds[jB] % jlev %
                                diffA % diffB % pair.distanceInKm %
                                corr % std::exp(expArg) % pges(jA, jlev) % pges(jB, jlev) % z;
        }  // verbose
      }  // jlev
    }  // single length scale or 2 length scales
  }  // pairs
}


void MetOfficeBuddyCheck::checkVectorData(const std::vector<MetOfficeBuddyPair> &pairs,
                                          const std::vector<int> &flags,
                                          const std::vector<bool> &verbose,
                                          const std::vector<float> &bgErrorHorizCorrScales,
                                          const std::vector<int> &stationIds,
                                          const std::vector<util::DateTime> &datetimes,
                                          const Eigen::ArrayXXf *pressures,
                                          const Eigen::ArrayXXf &uObsValues,
                                          const Eigen::ArrayXXf &vObsValues,
                                          const Eigen::ArrayXXf &obsErrors,
                                          const Eigen::ArrayXXf &uBgValues,
                                          const Eigen::ArrayXXf &vBgValues,
                                          const Eigen::ArrayXXf &bgErrors,
                                          Eigen::ArrayXXf &pges) const {
  using util::sqr;
  const boost::optional<int> &nolevs = options_.numLevels.value();
  const int numLevels = nolevs ? *nolevs : 0;  // number of actual levels
  const Eigen::Index cols = uObsValues.cols();  // number of data columns (levels) to loop over

  const bool isMaster = obsdb_.comm().rank() == 0;
  if (isMaster) {
    oops::Log::trace() << "ObsA  ObsB  StatIdA  StatIdB  lev  LDiffA LDiffB TDiffA TDiffB "
                          "Dist   Corr  Agree   PgeA   PgeB   Mult\n";
  }
  const double invTemporalCorrScale = 1.0 / options_.temporalCorrelationScale.value().toSeconds();
  const float missing = util::missingValue<float>();

  // Loop over buddy pairs
  for (const MetOfficeBuddyPair &pair : pairs) {
    const size_t jA = pair.obsIdA;
    const size_t jB = pair.obsIdB;

    if (!options_.useAllObservations) {
      // Check that observations are valid and buddy check is required
      if (flags[jA] != QCflags::pass || flags[jB] != QCflags::pass)
        continue;  // skip to next pair
    }

    // eqn 3.9
    // - hcScale: horizontal error scale for the pair of obs
    // - scaledDist: scaled Distance between observations
    const double hcScale = 0.5 * (bgErrorHorizCorrScales[jA] + bgErrorHorizCorrScales[jB]);
    const double scaledDist = pair.distanceInKm / hcScale;

    // Background error correlation between observation positions.
    // eqns 3.10, 3.11
    double lCorr;
    if (numLevels == 1) {
      lCorr = std::exp(-scaledDist -
                       options_.verticalCorrelationScale.value() *
                       sqr(std::log((*pressures)(jA, 0)/(*pressures)(jB, 0))) -
                       sqr((datetimes[jA] - datetimes[jB]).toSeconds() * invTemporalCorrScale));
    } else {
      // Multi-level/surface data; treat vertical correlation as 1.0
      lCorr = std::exp(-scaledDist -
                       sqr((datetimes[jA] - datetimes[jB]).toSeconds() * invTemporalCorrScale));
    }
    if ((1.0 + scaledDist) * lCorr < 0.1)  // Check against minimum background error correlation.
      continue;  // skip to next pair

    // Calculate longitudinal and transverse wind components
    const double sinRotA = std::sin(pair.rotationAInRad);
    const double cosRotA = std::cos(pair.rotationAInRad);
    const double sinRotB = std::sin(pair.rotationBInRad);
    const double cosRotB = std::cos(pair.rotationBInRad);

    // Loop over each level
    for (Eigen::Index jlev=0; jlev < cols; jlev++) {
      if (pges(jA, jlev) >= maxGrossErrorProbability ||
          pges(jB, jlev) >= maxGrossErrorProbability ||
          pges(jA, jlev) == missing || pges(jB, jlev) == missing)
        continue;  // skip to next level

      // Difference from background - longitudinal wind (eqn 3.19)
      const double lDiffA = cosRotA  * (uObsValues(jA, jlev) - uBgValues(jA, jlev))
          + sinRotA * (vObsValues(jA, jlev) - vBgValues(jA, jlev));
      // Difference from background - transverse wind (eqn 3.20)
      const double tDiffA = - sinRotA * (uObsValues(jA, jlev) - uBgValues(jA, jlev))
          + cosRotA * (vObsValues(jA, jlev) - vBgValues(jA, jlev));
      // Difference from background - longitudinal wind (eqn 3.19)
      const double lDiffB = cosRotB  * (uObsValues(jB, jlev) - uBgValues(jB, jlev))
          + sinRotB * (vObsValues(jB, jlev) - vBgValues(jB, jlev));
      // Difference from background - transverse wind (eqn 3.20)
      const double tDiffB = - sinRotB * (uObsValues(jB, jlev) - uBgValues(jB, jlev))
          + cosRotB  * (vObsValues(jB, jlev) - vBgValues(jB, jlev));

      // Estimated error variances (ob + bk; component wind variance; eqn 2.5)
      const double errVarA = sqr(obsErrors(jA, jlev)) + sqr(bgErrors(jA, jlev));
      const double errVarB = sqr(obsErrors(jB, jlev)) + sqr(bgErrors(jB, jlev));

      // Calculate covariances and probabilities eqn 3.12, 3.13
      const double lCovar = lCorr * bgErrors(jA, jlev) * bgErrors(jB, jlev);
      const double tCovar = (1.0 - options_.nonDivergenceConstraint * scaledDist) * lCovar;
      // rho2 = (total error correlation between ob positions)**2
      const double lRho2 = sqr(lCovar) / (errVarA * errVarB);  // eqn 3.14
      const double tRho2 = sqr(tCovar) / (errVarA * errVarB);  // eqn 3.14
      // Argument for exponents
      double expArg;
      if (std::abs (tRho2) <= 0.00001)
        expArg = 0.0;    // prevent division by tCovar=0.0
      else
        expArg = -(0.5 * tRho2 / (1.0 - tRho2)) *
            (sqr(tDiffA) / errVarA + sqr(tDiffB) / errVarB - 2.0 * tDiffA * tDiffB / tCovar);
      expArg = expArg - (0.5 * lRho2 / (1.0 - lRho2)) *
          (sqr(lDiffA) / errVarA + sqr(lDiffB) / errVarB - 2.0 * lDiffA * lDiffB / lCovar);
      expArg = options_.dampingFactor1 * (-0.5 * std::log((1.0 - lRho2) * (1.0 - lRho2)) + expArg);
      expArg = std::min(expArgMax, std::max(-expArgMax, expArg));           // eqn 3.22
      // Z = P(OA)*P(OB)/P(OA and OB)
      double z = 1.0 / (1.0 - (1.0 - pges(jA, jlev)) * (1.0 - pges(jB, jlev)) *
          (1.0 - std::exp(expArg)));
      if (z <= 0.0)
        z = 1.0;  // rounding error control
      z = std::pow(z, options_.dampingFactor2);         // eqn 3.16
      pges(jA, jlev) *= z;                                     // eqn 3.17
      pges(jB, jlev) *= z;                                     // eqn 3.17

      if (isMaster && (verbose[jA] || verbose[jB])) {
        oops::Log::trace() << boost::format("%5d %5d %8d %8d %5d "
                                            "%6.1f %6.1f %6.1f %6.1f %6.1f "
                                            "%5.3f %6.3f %6.3f %6.3f %6.3f\n") %
                              jA % jB % stationIds[jA] % stationIds[jB] % jlev %
                              lDiffA % lDiffB % tDiffA % tDiffB % pair.distanceInKm %
                              lCorr % std::exp(expArg) % pges(jA, jlev) % pges(jB, jlev) % z;
      }
    }
  }
}

std::vector<size_t> MetOfficeBuddyCheck::getValidObservationIds(
    const std::vector<bool> & apply,
    const boost::optional<Eigen::ArrayXXi> & profileIndex) const {
  std::vector<bool> isValid = apply;
  unselectRejectedLocations(isValid, filtervars_, *flags_,
                            UnselectLocationIf::ALL_FILTER_VARIABLES_REJECTED);

  std::vector<int> isValidAsInt(apply.begin(), apply.end());
  obsdb_.distribution()->allGatherv(isValidAsInt);
  isValid.assign(isValidAsInt.begin(), isValidAsInt.end());

  std::vector<size_t> validObsIds;

  const std::vector<int> isLevel1Valid = extract1stLev(isValidAsInt, profileIndex);
  for (size_t obsId = 0; obsId < isLevel1Valid.size(); ++obsId)
    if (isLevel1Valid[obsId])
      validObsIds.push_back(obsId);
  return validObsIds;
}

void MetOfficeBuddyCheck::flagRejectedObservations(
    const Variables &filtervars,
    const std::map<std::string, std::vector<float>> &grossErrProbsByVarName,
    std::vector<std::vector<bool>> &flagged) const {
  ASSERT(filtervars.size() == flagged.size());

  for (size_t varIndex = 0; varIndex < filtervars.size(); ++varIndex) {
    const std::string varName = fullVariableName(filtervars[varIndex]);
    const std::vector<float> &grossErrProbs = grossErrProbsByVarName.at(varName);
    std::vector<bool> &variableFlagged = flagged[varIndex];
    ASSERT(grossErrProbs.size() == variableFlagged.size());

    for (size_t obsId = 0; obsId < grossErrProbs.size(); ++obsId)
      if (grossErrProbs[obsId] >= options_.rejectionThreshold)
        variableFlagged[obsId] = true;
  }
}

void MetOfficeBuddyCheck::print(std::ostream & os) const {
  os << "MetOfficeBuddyCheck: config = " << options_ << std::endl;
}

}  // namespace ufo
