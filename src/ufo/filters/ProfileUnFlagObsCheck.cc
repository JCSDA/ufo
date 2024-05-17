/*
 * (C) British Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ProfileUnFlagObsCheck.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

/// \brief Helper classes to manage the optional piecewise linear scaling of the tolerance.

class AbsoluteTolerance {
 public:
  virtual ~AbsoluteTolerance() {}
  virtual float operator() (int job) const = 0;
};

class AbsoluteToleranceConstant : public AbsoluteTolerance {
 public:
  explicit AbsoluteToleranceConstant(const float tol) : unscaledTol_(tol) {}
  float operator() (int job) const override { return unscaledTol_; }
  static std::string classname() {return "AbsoluteToleranceConstant";}

 private:
  float unscaledTol_;
};

class AbsoluteTolerancePiecewise : public AbsoluteTolerance {
 public:
  AbsoluteTolerancePiecewise(const float tol, const std::map<float, float>& interpolationPoints,
                    const std::vector<float>& vertical) :
     unscaledTol_(tol), verticalData_(vertical) {
    std::vector<double> abscissas, ordinates;
    for (const std::pair<const float, float> &xy : interpolationPoints) {
      abscissas.push_back(xy.first);
      ordinates.push_back(xy.second);
    }
    toleranceScaleInterp = std::unique_ptr<PiecewiseLinearInterpolation>( new
                              PiecewiseLinearInterpolation{std::move(abscissas),
                                                           std::move(ordinates)});
  }

  float operator() (int job) const override {
      return unscaledTol_ * (*toleranceScaleInterp)(verticalData_[job]);
  }
  static std::string classname() {return "AbsoluteTolerancePiecewise";}

 private:
  const float unscaledTol_;
  const std::vector<float> verticalData_;
  std::unique_ptr<PiecewiseLinearInterpolation> toleranceScaleInterp;
};

class AbsoluteToleranceCreator {
 public:
  static void registerType(std::string name, AbsoluteToleranceCreator *factory) {
    getFactory()[name] = factory;
  }
  virtual std::unique_ptr<AbsoluteTolerance> create_unique(
    const ProfileUnFlagObsCheck::Parameters_& parameters, const ioda::ObsSpace& obsdb) = 0;
  static std::unique_ptr<AbsoluteTolerance> create_unique(std::string name,
    const ProfileUnFlagObsCheck::Parameters_& parameters, const ioda::ObsSpace& obsdb) {
    std::unique_ptr<AbsoluteTolerance> absoluteTolerance =
        std::move(getFactory()[name]->create_unique(parameters, obsdb));
    return absoluteTolerance;
  }
  static std::map<std::string, AbsoluteToleranceCreator*> &getFactory() {
    static std::map<std::string, AbsoluteToleranceCreator*> creatorMap;
    return creatorMap;
  }
};

class AbsoluteToleranceConstantCreator : public AbsoluteToleranceCreator {
 public:
  AbsoluteToleranceConstantCreator()
    { AbsoluteToleranceCreator::registerType(AbsoluteToleranceConstant::classname(), this); }
  std::unique_ptr<AbsoluteTolerance> create_unique(
      const ProfileUnFlagObsCheck::Parameters_& parameters, const ioda::ObsSpace& obsdb) {
    std::unique_ptr<AbsoluteTolerance> absTolConst(
        new AbsoluteToleranceConstant(parameters.aTol.value()));
    return absTolConst;
  }
};
static AbsoluteToleranceConstantCreator absoluteToleranceConstantCreator;

class AbsoluteTolerancePiecewiseCreator : public AbsoluteToleranceCreator {
 public:
  AbsoluteTolerancePiecewiseCreator()
    { AbsoluteToleranceCreator::registerType(AbsoluteTolerancePiecewise::classname(), this); }
  std::unique_ptr<AbsoluteTolerance> create_unique(
      const ProfileUnFlagObsCheck::Parameters_& parameters, const ioda::ObsSpace& obsdb) {
    std::vector<float> verticalData_(obsdb.nlocs());
    const Variable verticalCoord = parameters.verticalCoord.value();
    obsdb.get_db(verticalCoord.group(), verticalCoord.variable(), verticalData_);

    std::unique_ptr<AbsoluteTolerance> absTolPiecewise(
        new AbsoluteTolerancePiecewise(
              parameters.aTol.value(),
              parameters.verticalToleranceScaleInterpolationPoints.value().value(),
              verticalData_));
    return absTolPiecewise;
  }
};
static AbsoluteTolerancePiecewiseCreator absoluteTolerancePiecewiseCreator;

// -----------------------------------------------------------------------------

ProfileUnFlagObsCheck::ProfileUnFlagObsCheck(
        ioda::ObsSpace & obsdb,
        const Parameters_ & parameters,
        std::shared_ptr<ioda::ObsDataVector<int> > flags,
        std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "ProfileUnFlagObsCheck constructor" << std::endl;
}

// -----------------------------------------------------------------------------

ProfileUnFlagObsCheck::~ProfileUnFlagObsCheck() {
  oops::Log::trace() << "ProfileUnFlagObsCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/// Apply the profile check for the number of observations.

void ProfileUnFlagObsCheck::applyFilter(const std::vector<bool> & apply,
                                        const Variables & filtervars,
                                        std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ProfileUnFlagObsCheck preProcess filter" << std::endl;
  const oops::ObsVariables observed = obsdb_.obsvariables();

  // setup the calculation of the absolute tolerance to use a piecewise linear
  // function if specified.
  std::unique_ptr<AbsoluteTolerance> absoluteToleranceCalculate;

  if (parameters_.verticalToleranceScaleInterpolationPoints.value()) {
    absoluteToleranceCalculate = std::move(AbsoluteToleranceCreator::create_unique(
          AbsoluteTolerancePiecewise::classname(), parameters_, obsdb_));
  } else {
    absoluteToleranceCalculate = std::move(AbsoluteToleranceCreator::create_unique(
          AbsoluteToleranceConstant::classname(), parameters_, obsdb_));
  }

  // Check the number of channels and variables to process
  const size_t nChans = std::max(obsdb_.nchans(), 1LU);
  // For multi-level data "nvars" is the number of simulated variables times the
  // number of channels.  Therefore, divide "nvars" by the number of channels
  // to get the number of simulated variables.
  const size_t nActualVars = filtervars.nvars() / nChans;

  std::vector<float> variableData(obsdb_.nlocs());

  // Get the record numbers from the observation data.  These will be used to identify
  // which observations belong to which profile.
  const std::vector<size_t> & recordNumbers = obsdb_.recidx_all_recnums();

  // Loop over the number of actual variables
  for (size_t iVar = 0; iVar < nActualVars; ++iVar) {
    // Loop over the unique profiles to check for valid data and unflag if there is agreement
    obsdb_.get_db(filtervars.variable(iVar).group(), filtervars.variable(iVar).variable(),
                  variableData);

    std::vector<size_t> variableIndicesMap(nChans, -1);
    for (size_t iChan = 0; iChan < nChans; ++iChan) {
        const size_t iFilterVar = iVar * nChans + iChan;
        const std::string variableName = filtervars.variable(iFilterVar).variable();
        variableIndicesMap[iChan] = observed.find(variableName);
    }

    for (size_t iProfile : recordNumbers) {
      const std::vector<size_t> & allObsNumbers = obsdb_.recidx_vector(iProfile);

      // Retrieve indices of profiles/observations where we intend to apply unflagging
      std::vector<size_t> obsNumbers;
      for (size_t jObs : allObsNumbers)
        if (apply[jObs])
          obsNumbers.push_back(jObs);

      // Accept observations with valid neighbours within a tolerance.
      for (size_t iChan = 0; iChan < nChans; ++iChan) {
        const size_t iFilterVar = iVar * nChans + iChan;
        const size_t jVar = variableIndicesMap[iChan];
        int numMadeValid = 0;
        for (size_t jObsProfile = 0; jObsProfile < obsNumbers.size(); ++jObsProfile) {
          size_t jObsGlobal = obsNumbers[jObsProfile];
          const float absTol = (*absoluteToleranceCalculate)(jObsGlobal);
          if (((*flags_)[jVar][jObsGlobal] != QCflags::pass) &&
              ((*flags_)[jVar][jObsGlobal] != QCflags::passive) &&
              ((*flags_)[jVar][jObsGlobal] != QCflags::missing)) {
            bool makeValid = false;
            if (jObsProfile == 0) {
              makeValid = true;
            } else {
              const size_t jObsMinus = obsNumbers[jObsProfile-1];
              makeValid = (*flags_)[jVar][jObsMinus] == QCflags::pass &&
                  (std::abs(variableData[jObsMinus] - variableData[jObsGlobal]) < absTol);
            }
            if (jObsProfile < obsNumbers.size()-1) {
              const size_t jObsPlus = obsNumbers[jObsProfile+1];
              makeValid = makeValid &&
                  ((*flags_)[jVar][jObsPlus] == QCflags::pass) &&
                  (std::abs(variableData[jObsPlus] - variableData[jObsGlobal]) < absTol);
            }
            if (makeValid) {
              numMadeValid++;
              flagged[iFilterVar][jObsGlobal] = true;
            }
          }
        }
        if ((iChan == 0) && (numMadeValid > 0))
          oops::Log::debug() << "\nUnFlagObsCheck " << observed[variableIndicesMap[iChan]]
            << " profile "  << iProfile << " [channel] numMadeValid:";
        if (numMadeValid > 0)
          oops::Log::debug() << " [" << iChan << "] " << numMadeValid;
      }
    }
  }

  oops::Log::debug() << std::endl;
}

// -----------------------------------------------------------------------------

void ProfileUnFlagObsCheck::print(std::ostream & os) const {
  os << "ProfileUnFlagObsCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
