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
#include <string>
#include <utility>
#include <vector>

#include <boost/format.hpp>

#include "eckit/config/Configuration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/sqr.h"
#include "ufo/filters/MetOfficeBuddyCheckParameters.h"
#include "ufo/filters/MetOfficeBuddyPair.h"
#include "ufo/filters/MetOfficeBuddyPairFinder.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

namespace {

const double expArgMax = 80.0;  // maximum argument of the exponential function
const float maxGrossErrorProbability = 1.0;  // maximum allowed value for PGE

std::string fullVariableName(const Variable &var)
{
  if (var.group().empty())
    return var.variable();
  else
    return var.variable() + "@" + var.group();
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

struct ScalarSingleLevelVariableData {
  const std::vector<int> *varFlags = nullptr;
  const std::vector<float> *obsValues = nullptr;
  const std::vector<float> *obsBiases = nullptr;
  const std::vector<float> *obsErrors = nullptr;
  std::vector<float> bgValues;
  std::vector<float> bgErrors;
  std::vector<float> grossErrorProbabilities;
};

}  // namespace

/// \brief Metadata of all observations processed by the filter.
struct MetOfficeBuddyCheck::MetaData {
  std::vector<float> latitudes;
  std::vector<float> longitudes;
  std::vector<util::DateTime> datetimes;
  boost::optional<std::vector<float>> pressures;
  std::vector<int> stationIds;
};

MetOfficeBuddyCheck::MetOfficeBuddyCheck(ioda::ObsSpace& obsdb, const Parameters_& parameters,
                                         std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                         std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, std::move(flags), std::move(obserr)), options_(parameters)
{
  oops::Log::debug() << "MetOfficeBuddyCheck: config = " << options_ << std::endl;
  allvars_ += Variables(filtervars_, "HofX");
}

// Required for the correct destruction of options_.
MetOfficeBuddyCheck::~MetOfficeBuddyCheck()
{}

void MetOfficeBuddyCheck::applyFilter(const std::vector<bool> & apply,
                                      const Variables & filtervars,
                                      std::vector<std::vector<bool>> & flagged) const {
  // Identify observations to process and collect their data and metadata

  const std::vector<size_t> validObsIds = getValidObservationIds(apply);
  MetaData obsData = collectMetaData();
  const std::vector<float> bgErrorHorizCorrScales = calcBackgroundErrorHorizontalCorrelationScales(
        validObsIds, obsData.latitudes);
  const std::vector<bool> verbose = flagAndPrintVerboseObservations(
        validObsIds, obsData.latitudes, obsData.longitudes, obsData.datetimes,
        obsData.pressures.get_ptr(), obsData.stationIds, bgErrorHorizCorrScales);

  const oops::Variables &observedVars = obsdb_.obsvariables();
  ioda::ObsDataVector<float> obsValues(obsdb_, filtervars.toOopsVariables(), "ObsValue");
  ioda::ObsDataVector<float> obsBiases(obsdb_, filtervars.toOopsVariables(), "ObsBias",
                                       false /*fail if not found?*/);

  auto getFilterVariableName = [&] (size_t filterVarIndex) {
    return filtervars.variable(filterVarIndex).variable();
  };

  auto getScalarSingleLevelVariableData = [&] (size_t filterVarIndex) {
    const size_t observedVarIndex = observedVars.find(getFilterVariableName(filterVarIndex));

    ScalarSingleLevelVariableData data;
    data.varFlags = &(*flags_)[observedVarIndex];
    data.obsValues = &obsValues[filterVarIndex];
    data.obsBiases = &obsBiases[filterVarIndex];
    data.obsErrors = &(*obserr_)[observedVarIndex];
    data_.get(ufo::Variable(filtervars[filterVarIndex], "HofX"), data.bgValues);
    // TODO(wsmigaj): It isn't clear where to get background error variances from.
    data_.get(ufo::Variable(filtervars[filterVarIndex], "HofXError"), data.bgErrors);
    // TODO(wsmigaj): How is this variable going to be initialized?
    data_.get(ufo::Variable(filtervars[filterVarIndex], "GrossErrorProbability"),
              data.grossErrorProbabilities);
    return data;
  };

  // Identify buddy pairs

  const std::vector<float> *pressures =
      options_.sortByPressure ? obsData.pressures.get_ptr() : nullptr;
  MetOfficeBuddyPairFinder buddyPairFinder(options_, obsData.latitudes, obsData.longitudes,
                                           obsData.datetimes, pressures,
                                           obsData.stationIds);
  const std::vector<MetOfficeBuddyPair> buddyPairs = buddyPairFinder.findBuddyPairs(validObsIds);

  // Buddy-check all filter variables

  // Gross error probabilities updated by buddy check, indexed by variable name.
  // (The updated values are copied to the ObsSpace only after all variables have been processed).
  std::map<std::string, std::vector<float>> calculatedGrossErrProbsByVarName;

  bool previousVariableWasFirstComponentOfTwo = false;
  for (size_t filterVarIndex = 0; filterVarIndex < filtervars.size(); ++filterVarIndex) {
    if (previousVariableWasFirstComponentOfTwo) {
      // Vector (two-component) variable
      ScalarSingleLevelVariableData firstComponentData =
          getScalarSingleLevelVariableData(filterVarIndex - 1);
      ScalarSingleLevelVariableData secondComponentData =
          getScalarSingleLevelVariableData(filterVarIndex);

      for (size_t i = 0; i < firstComponentData.grossErrorProbabilities.size(); ++i)
        firstComponentData.grossErrorProbabilities[i] =
            std::max(firstComponentData.grossErrorProbabilities[i],
                     secondComponentData.grossErrorProbabilities[i]);

      checkVectorSurfaceData(buddyPairs,
                             *firstComponentData.varFlags,
                             verbose, bgErrorHorizCorrScales,
                             obsData.stationIds, obsData.datetimes,
                             *firstComponentData.obsValues, *firstComponentData.obsBiases,
                             *secondComponentData.obsValues, *secondComponentData.obsBiases,
                             *firstComponentData.obsErrors,
                             firstComponentData.bgValues, secondComponentData.bgValues,
                             firstComponentData.bgErrors,
                             firstComponentData.grossErrorProbabilities);
      // OPS doesn't update the gross error probabilities of the second component variable,
      // but it seems more consistent to do so (and it facilitates the implementation of
      // flagRejectedObservations().
      // The implementation of checkVectorSurfaceData() still assumes that the *input* gross error
      // probabilities, flags and background error estimates are the same for both components.
      calculatedGrossErrProbsByVarName[getFilterVariableName(filterVarIndex - 1)] =
          firstComponentData.grossErrorProbabilities;
      calculatedGrossErrProbsByVarName[getFilterVariableName(filterVarIndex)] =
          std::move(firstComponentData.grossErrorProbabilities);

      previousVariableWasFirstComponentOfTwo = false;
    } else {
      if (filtervars[filterVarIndex].options().getBool("first_component_of_two", false)) {
        previousVariableWasFirstComponentOfTwo = true;
      } else {
        // Scalar variable
        ScalarSingleLevelVariableData data =
            getScalarSingleLevelVariableData(filterVarIndex);
        checkScalarSurfaceData(buddyPairs,
                               *data.varFlags, verbose, bgErrorHorizCorrScales,
                               obsData.stationIds, obsData.datetimes,
                               *data.obsValues, *data.obsBiases, *data.obsErrors,
                               data.bgValues, data.bgErrors,
                               data.grossErrorProbabilities);
        calculatedGrossErrProbsByVarName[getFilterVariableName(filterVarIndex)] =
            std::move(data.grossErrorProbabilities);
      }
    }
  }

  // Update observations flags and gross error probabilities

  for (const auto &varNameAndGrossErrProbs : calculatedGrossErrProbsByVarName) {
    obsdb_.put_db("GrossErrorProbability", varNameAndGrossErrProbs.first,
                  varNameAndGrossErrProbs.second);
  }

  flagRejectedObservations(filtervars, calculatedGrossErrProbsByVarName, flagged);

  if (filtervars.size() != 0) {
    oops::Log::trace() << "MetOfficeBuddyCheck: flagged? = " << flagged[0] << std::endl;
  }
}

MetOfficeBuddyCheck::MetaData MetOfficeBuddyCheck::collectMetaData() const {
  MetaData obsData;

  obsData.latitudes.resize(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "latitude", obsData.latitudes);

  obsData.longitudes.resize(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "longitude", obsData.longitudes);

  obsData.datetimes.resize(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "datetime", obsData.datetimes);

  if (obsdb_.has("MetaData", "air_pressure")) {
    obsData.pressures = std::vector<float>(obsdb_.nlocs());
    obsdb_.get_db("MetaData", "air_pressure", *obsData.pressures);
  }

  obsData.stationIds = getStationIds();

  return obsData;
}

std::vector<int> MetOfficeBuddyCheck::getStationIds() const {
  const boost::optional<Variable> &stationIdVariable = options_.stationIdVariable.value();
  if (stationIdVariable == boost::none) {
    if (obsdb_.obs_group_var().empty()) {
      // Observations were not grouped into records.
      // Assume all observations were taken by the same station.
      return std::vector<int>(obsdb_.nlocs(), 0);
    } else {
      const std::vector<size_t> &recordNumbers = obsdb_.recnum();
      return std::vector<int>(recordNumbers.begin(), recordNumbers.end());
    }
  } else {
    switch (obsdb_.dtype(stationIdVariable->group(), stationIdVariable->variable())) {
    case ioda::ObsDtype::Integer:
      {
        std::vector<int> stationIds(obsdb_.nlocs());
        obsdb_.get_db(stationIdVariable->group(), stationIdVariable->variable(), stationIds);
        return stationIds;
      }

    case ioda::ObsDtype::String:
      {
        std::vector<std::string> stringIds(obsdb_.nlocs());
        obsdb_.get_db(stationIdVariable->group(), stationIdVariable->variable(), stringIds);
        return mapDistinctValuesToDistinctInts(stringIds);
      }

    default:
      throw eckit::UserError("Only integer and string variables may be used as station IDs",
                             Here());
    }
  }
}

std::vector<float> MetOfficeBuddyCheck::calcBackgroundErrorHorizontalCorrelationScales(
    const std::vector<size_t> &validObsIds, const std::vector<float> &latitudes) const {

  std::vector<double> abscissas, ordinates;
  for (const std::pair<const float, float> &xy :
       options_.horizontalCorrelationScaleInterpolationPoints.value()) {
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

void MetOfficeBuddyCheck::checkScalarSurfaceData(const std::vector<MetOfficeBuddyPair> &pairs,
                                                 const std::vector<int> &flags,
                                                 const std::vector<bool> &verbose,
                                                 const std::vector<float> &bgErrorHorizCorrScales,
                                                 const std::vector<int> &stationIds,
                                                 const std::vector<util::DateTime> &datetimes,
                                                 const std::vector<float> &obsValues,
                                                 const std::vector<float> &obsBiases,
                                                 const std::vector<float> &obsErrors,
                                                 const std::vector<float> &bgValues,
                                                 const std::vector<float> &bgErrors,
                                                 std::vector<float> &pges) const {
  using util::sqr;

  const bool isMaster = obsdb_.comm().rank() == 0;
  if (isMaster) {
    oops::Log::trace() << __func__ << " "
                       << " dampingFactor1 = " << options_.dampingFactor1
                       << ", dampingFactor2 = " << options_.dampingFactor2 << '\n';
    oops::Log::trace() << "ObsA  ObsB  StatIdA  StatIdB  DiffA DiffB "
                          "Dist   Corr  Agree   PgeA   PgeB   Mult\n";
  }

  const double invTemporalCorrScale = 1.0 / options_.temporalCorrelationScale.value().toSeconds();

  for (const MetOfficeBuddyPair &pair : pairs) {
    const size_t jA = pair.obsIdA;
    const size_t jB = pair.obsIdB;

    // Check that observations are valid and buddy check is required
    if (!(flags[jA] == QCflags::pass && flags[jB] == QCflags::pass &&
          pges[jA] < maxGrossErrorProbability && pges[jB] < maxGrossErrorProbability))
      continue;

    // eqn 3.9
    const double hcScale = 0.5 * (bgErrorHorizCorrScales[jA] + bgErrorHorizCorrScales[jB]);
    const double scaledDist = pair.distanceInKm / hcScale;

    // Background error correlation between ob positions.
    // Surface data; treat vertical correlation as 1.0
    // eqns 3.10, 3.11
    const double corr = (1.0 + scaledDist) *
        std::exp(-scaledDist - sqr((datetimes[jA] - datetimes[jB]).toSeconds() *
                                   invTemporalCorrScale));

    if (corr < 0.1)
      continue;  // skip to next pair

    // Differences from background
    double diffA = obsValues[jA] + obsBiases[jA] - bgValues[jA];
    double diffB = obsValues[jB] + obsBiases[jB] - bgValues[jB];
    // Estimated error variances (ob+bk) (eqn 2.5)
    double errVarA = sqr(obsErrors[jA]) + sqr(bgErrors[jA]);
    double errVarB = sqr(obsErrors[jB]) + sqr(bgErrors[jB]);
    // Background error covariance between ob positions (eqn 3.13)
    double covar = corr * bgErrors[jA] * bgErrors[jB];
    // (Total error correlation between ob positions)**2 (eqn 3.14)
    double rho2 = sqr(covar) / (errVarA * errVarB);
    // Argument for exponents
    double expArg = -(0.5 * rho2 / (1.0 - rho2)) *
        (sqr(diffA) / errVarA + sqr(diffB) / errVarB - 2.0 * diffA * diffB / covar);
    expArg = options_.dampingFactor1 * (-0.5 * std::log(1.0 - rho2) + expArg);  // exponent of
    expArg = std::min(expArgMax, std::max(-expArgMax, expArg));                 // eqn 3.18
    // Z = P(OA)*P(OB)/P(OA and OB)
    double z = 1.0 / (1.0 - (1.0 - pges[jA]) * (1.0 - pges[jB]) * (1.0 - std::exp(expArg)));
    if (z <= 0.0)
      z = 1.0;  // rounding error control
    z = std::pow(z, options_.dampingFactor2);  // eqn 3.16
    pges[jA] *= z;                              // eqn 3.17
    pges[jB] *= z;                              // eqn 3.17
    if (isMaster && (verbose[jA] || verbose[jB])) {
      oops::Log::trace() << boost::format("%5d %5d %8d %8d "
                                          "%5.1f %5.1f %6.1f "
                                          "%5.3f %6.3f %6.3f %6.3f %6.3f\n") %
                            jA % jB % stationIds[jA] % stationIds[jB] %
                            diffA % diffB % pair.distanceInKm %
                            corr % std::exp(expArg) % pges[jA] % pges[jB] % z;
    }
  }
}

void MetOfficeBuddyCheck::checkVectorSurfaceData(const std::vector<MetOfficeBuddyPair> &pairs,
                                                 const std::vector<int> &flags,
                                                 const std::vector<bool> &verbose,
                                                 const std::vector<float> &bgErrorHorizCorrScales,
                                                 const std::vector<int> &stationIds,
                                                 const std::vector<util::DateTime> &datetimes,
                                                 const std::vector<float> &uObsValues,
                                                 const std::vector<float> &uObsBiases,
                                                 const std::vector<float> &vObsValues,
                                                 const std::vector<float> &vObsBiases,
                                                 const std::vector<float> &obsErrors,
                                                 const std::vector<float> &uBgValues,
                                                 const std::vector<float> &vBgValues,
                                                 const std::vector<float> &bgErrors,
                                                 std::vector<float> &pges) const {
  using util::sqr;

  const bool isMaster = obsdb_.comm().rank() == 0;
  if (isMaster) {
    oops::Log::trace() << __func__ << " "
                       << " dampingFactor1 = " << options_.dampingFactor1
                       << ", dampingFactor2 = " << options_.dampingFactor2 << '\n';
    oops::Log::trace() << "ObsA  ObsB  StatIdA  StatIdB  LDiffA LDiffB TDiffA TDiffB "
                          "Dist   Corr  Agree   PgeA   PgeB   Mult\n";
  }

  const double invTemporalCorrScale = 1.0 / options_.temporalCorrelationScale.value().toSeconds();

  for (const MetOfficeBuddyPair &pair : pairs) {
    const size_t jA = pair.obsIdA;
    const size_t jB = pair.obsIdB;

    // Check that observations are valid and buddy check is required
    if (!(flags[jA] == QCflags::pass && flags[jB] == QCflags::pass &&
          pges[jA] < maxGrossErrorProbability && pges[jB] < maxGrossErrorProbability))
      continue;

    // eqn 3.9
    double horizCorrScale = 0.5 * (bgErrorHorizCorrScales[jA] + bgErrorHorizCorrScales[jB]);
    double scaleDist = pair.distanceInKm / horizCorrScale;
    // Background error correlation between ob positions.
    // Surface data; treat vertical correlation as 1.0
    // eqns 3.10, 3.11
    const double lCorr = std::exp(-scaleDist - sqr((datetimes[jA] - datetimes[jB]).toSeconds() *
                                   invTemporalCorrScale));

    if ((1.0 + scaleDist) * lCorr < 0.1)
      continue;  // skip to next pair

    // Calculate longitudinal and transverse wind components
    double sinRot = std::sin(pair.rotationAInRad);
    double cosRot = std::cos(pair.rotationAInRad);
    // Difference from background - longitudinal wind
    double lDiffA = cosRot  * (uObsValues[jA] + uObsBiases[jA] - uBgValues[jA])
        + sinRot * (vObsValues[jA] + vObsBiases[jA] - vBgValues[jA]);           // eqn 3.19
    // Difference from background - transverse wind
    double tDiffA = - sinRot * (uObsValues[jA] + uObsBiases[jA] - uBgValues[jA])
        + cosRot * (vObsValues[jA] + vObsBiases[jA] - vBgValues[jA]);           // eqn 3.20
    sinRot = std::sin(pair.rotationBInRad);
    cosRot  = std::cos(pair.rotationBInRad);
    // Difference from background - longitudinal wind
    double lDiffB = cosRot  * (uObsValues[jB] + uObsBiases[jB] - uBgValues[jB])
        + sinRot * (vObsValues[jB] + vObsBiases[jB] - vBgValues[jB]);           // eqn 3.19
    // Difference from background - transverse wind
    double tDiffB = - sinRot * (uObsValues[jB] + uObsBiases[jB] - uBgValues[jB])
        + cosRot  * (vObsValues[jB] + vObsBiases[jB] - vBgValues[jB]);          // eqn 3.20

    // Estimated error variances (ob + bk; component wind variance)
    double errVarA = sqr(obsErrors[jA]) + sqr(bgErrors[jA]);                // eqn 2.5
    double errVarB = sqr(obsErrors[jB]) + sqr(bgErrors[jB]);                // eqn 2.5

    // Calculate covariances and probabilities
    double lCovar = lCorr * bgErrors[jA] * bgErrors[jB];                    // eqn 3.13
    double tCovar = (1.0 - options_.nonDivergenceConstraint * scaleDist) * lCovar;  // eqn 3.12, 13
    // rho2 = (total error correlation between ob positions)**2
    double lRho2 = sqr(lCovar) / (errVarA * errVarB);                       // eqn 3.14
    double tRho2 = sqr(tCovar) / (errVarA * errVarB);                       // eqn 3.14
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
    double z = 1.0 / (1.0 - (1.0 - pges[jA]) * (1.0 - pges[jB]) * (1.0 - std::exp(expArg)));
    if (z <= 0.0)
      z = 1.0;  // rounding error control
    z = std::pow(z, options_.dampingFactor2);         // eqn 3.16
    pges[jA] *= z;                                     // eqn 3.17
    pges[jB] *= z;                                     // eqn 3.17

    if (isMaster && (verbose[jA] || verbose[jB])) {
      oops::Log::trace() << boost::format("%5d %5d %8d %8d "
                                          "%6.1f %6.1f %6.1f %6.1f %6.1f "
                                          "%5.3f %6.3f %6.3f %6.3f %6.3f\n") %
                            jA % jB % stationIds[jA] % stationIds[jB] %
                            lDiffA % lDiffB % tDiffA % tDiffB % pair.distanceInKm %
                            lCorr % std::exp(expArg) % pges[jA] % pges[jB] % z;
    }
  }
}

std::vector<size_t> MetOfficeBuddyCheck::getValidObservationIds(
    const std::vector<bool> & apply) const {
  std::vector<size_t> validObsIds;
  for (size_t obsId = 0; obsId < apply.size(); ++obsId)
    // TODO(wsmigaj): The second condition below may need reviewing.
    // Perhaps we should process only observations marked as passed in all filter variables?
    // Or those marked as passed in at least one filter variable?
    if (apply[obsId] && (*flags_)[0][obsId] == QCflags::pass)
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
