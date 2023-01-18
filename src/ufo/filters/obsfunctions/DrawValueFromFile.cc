/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>

#include "eckit/exception/Exceptions.h"
#include "ioda/Misc/StringFuncs.h"  // for convertV1PathToV2Path
#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/DrawValueFromFile.h"

namespace ufo {

namespace {

// -----------------------------------------------------------------------------
typedef std::vector<std::pair<std::string, boost::variant<std::vector<int>,
                                                          std::vector<float>,
                                                          std::vector<std::string>,
                                                          std::vector<std::vector<int>>,
                                                          std::vector<std::vector<float>>,
                                                          std::vector<std::vector<std::string>>>
                   >> ObData;


// -----------------------------------------------------------------------------
/// A boost::variant visitor whose function call operator prints the value of the specified element
/// of a vector to the specified output stream.
class PrintVisitor : public boost::static_visitor<void> {
 public:
  PrintVisitor(std::ostream &os, size_t iloc) :
    os_(os), iloc_(iloc) {}

  template <typename T>
  void operator()(const std::vector<T> &obData) {
    os_ << obData[iloc_];
  }

  template <typename T>
  void operator()(const std::vector<std::vector<T>> &obData) {
    for (size_t i=0; i < obData[iloc_].size(); ++i) {
      os_ << obData[iloc_][i] << "  ";
    }
  }

 private:
  std::ostream &os_;
  size_t iloc_;
};


// -----------------------------------------------------------------------------
/// \brief This is a convenience function for updating our container for useful observation data
template <typename T>
void updateObData(const ObsFilterData &in, const Variable &var, ObData &obData) {
  std::vector<T> dat;
  in.get(var, dat);
  obData.emplace_back(var.fullName(), std::move(dat));
}

// -----------------------------------------------------------------------------
/// \brief This is a convenience function for updating our container for useful
/// data, given as an input vector
template <typename T>
void updateObData(const Variable &var, const std::vector<T> &dat, ObData &obData) {
  obData.emplace_back(var.fullName(), std::move(dat));
}

// -----------------------------------------------------------------------------
/// \brief This is a convenience function for updating our container for useful observation data
template <typename T>
void updateObData(const ObsFilterData &in, const Variable &var, ObData &obData,
                  const std::vector<int> &channels) {
  std::vector<std::vector<T>> dat;
  for (size_t ichan=0; ichan < channels.size(); ++ichan) {
    std::vector<T> temp;
    in.get(Variable(var.fullName(), channels)[ichan], temp);
    dat.push_back(temp);
  }
  obData.emplace_back(var.fullName(), std::move(dat));
}

// -----------------------------------------------------------------------------
/// \brief Add datetime observation information data to our container.
/// \details We simply convert the datetimes to strings as our implementation does not discriminate
/// between the two types.
void updateObDataDateTime(const ObsFilterData &in, const Variable &var, ObData &obData) {
  std::vector<util::DateTime> dat;
  std::vector<std::string> datConv;
  in.get(var, dat);
  datConv.resize(dat.size());

  // Convert the vec. of datetime. to strings
  std::transform(dat.begin(), dat.end(), datConv.begin(),
                 [](util::DateTime dt){return dt.toString();});
  obData.emplace_back(var.fullName(), std::move(datConv));
}

}  // anonymous namespace

// -----------------------------------------------------------------------------
constexpr char InterpMethodParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<InterpMethod>
  InterpMethodParameterTraitsHelper::namedValues[];

// -----------------------------------------------------------------------------
constexpr char ExtrapolationModeParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<ExtrapolationMode>
  ExtrapolationModeParameterTraitsHelper::namedValues[];

// -----------------------------------------------------------------------------
constexpr char EquidistantChoiceParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<EquidistantChoice>
  EquidistantChoiceParameterTraitsHelper::namedValues[];

// -----------------------------------------------------------------------------
constexpr char CoordinateTransformationParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<CoordinateTransformation>
  CoordinateTransformationParameterTraitsHelper::namedValues[];

// -----------------------------------------------------------------------------
static ObsFunctionMaker<DrawValueFromFile<float>> floatMaker("DrawValueFromFile");
static ObsFunctionMaker<DrawValueFromFile<int>> intMaker("DrawValueFromFile");
static ObsFunctionMaker<DrawValueFromFile<std::string>> stringMaker("DrawValueFromFile");

// -----------------------------------------------------------------------------
// Could be moved to a non-templated base class
template <typename T>
DrawValueFromFile<T>::DrawValueFromFile(const eckit::LocalConfiguration &config)
  : allvars_() {
    // Initialize options
    options_.deserialize(config);

    std::vector<eckit::LocalConfiguration> interpSubConfs;
    const std::vector<InterpolationParameters> &interpolationParameters =
      options_.interpolation.value();

    int nlin = 0;
    int nbilin = 0;
    int ntrilin = 0;
    int nlintotal = 0;
    for (auto intParam = interpolationParameters.begin();
         intParam != interpolationParameters.end(); ++intParam) {
      const ufo::InterpMethod & method = intParam->method.value();
      if (method == InterpMethod::TRILINEAR) {
        ntrilin++;
        nlintotal++;
      } else if (method == InterpMethod::BILINEAR) {
        nbilin++;
        nlintotal++;
      } else if (method == InterpMethod::LINEAR) {
        nlin++;
        nlintotal++;
      } else if (ntrilin > 0) {
        // if trilinear interpolation has been used at least once
        // but the current option is not trilinear interpolation
        throw eckit::UserError("Trilinear interpolation can only be supplied as the final three "
                               "arguments.", Here());
      } else if (nbilin > 0) {
        // if bilinear interpolation has been used at least once
        // but the current option is not bilinear interpolation
        throw eckit::UserError("Bilinear interpolation can only be supplied as the final two "
                               "arguments.", Here());
      } else if (nlin > 0) {
        // if linear interpolation has been used at least once
        // but the current option is not linear interpolation
        throw eckit::UserError("Linear interpolation can only be supplied as the very last "
                               "argument.", Here());
      }
      // if more than one interpolation option has been used
      if ((ntrilin > 0 && ntrilin != nlintotal) ||
          (nbilin > 0 && nbilin != nlintotal) ||
          (nlin > 0 && nlin != nlintotal)) {
        throw eckit::UserError("Cannot use multiple interpolation methods.", Here());
      }

      const std::string varName = ioda::convertV1PathToV2Path(intParam->name.value());
      interpSubConfs.push_back(intParam->toConfiguration());
      interpMethod_[varName] = method;
      extrapMode_[varName] = intParam->extrapMode.value();
      equidistantChoice_[varName] = intParam->equidistanceChoice.value();
      useChannelList_[varName] = intParam->useChannelList.value();
      coordinateTransformation_[varName] = intParam->coordinateTransformation.value();
    }
    // if linear interpolation has been configured incorrectly
    if (nlin > 0 && nlin != 1) {
      throw eckit::UserError("Linear interpolation requires one variable.", Here());
    }
    // if bilinear interpolation has been configured incorrectly
    if (nbilin > 0 && nbilin != 2) {
      throw eckit::UserError("Bilinear interpolation requires two variables.", Here());
    }
    // if trilinear interpolation has been configured incorrectly
    if (ntrilin > 0 && ntrilin != 3) {
      throw eckit::UserError("Trilinear interpolation requires three variables.", Here());
    }

    // Get channels from options
    std::set<int> channelset = oops::parseIntSet(options_.chlist);
    std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));

    fpath_ = options_.fpath.value();
    allvars_  = Variables(interpSubConfs);
}


// -----------------------------------------------------------------------------
template <typename ExtractedValue>
class ExtractVisitor : public boost::static_visitor<void> {
 public:
  ExtractVisitor(DataExtractor<ExtractedValue> &interpolator, size_t iloc) :
    interpolator(interpolator), iloc(iloc), ivar(0) {}

  // Get the correct element for a vector
  template <typename T>
  T getElement(const std::vector<T> &vec) {
    return vec[iloc];
  }

  // Get the correct element for a vector of vectors
  template <typename T>
  T getElement(const std::vector<std::vector<T>> &vec) {
    return vec[ivar][iloc];
  }

  // Catch-all function for tri-linear interpolation
  template <typename T>
  float getFloatElement(const T &vec) {
    throw eckit::UserError("Trilinear interpolation must be performed with float coordinates.",
                           Here());
  }

  // Get the correct element for a float for use with trilinear interpolation
  float getFloatElement(const std::vector<float> &vec) {
    return getElement(vec);
  }

  // Get the correct element for a float for use with trilinear interpolation
  float getFloatElement(const std::vector<std::vector<float>> &vec) {
    return getElement(vec);
  }

  // Extract data for 1D extraction (nearest, linear interpolation etc.)
  template <typename T>
  void operator()(const T &obDat) {
    interpolator.extract(getElement(obDat));
  }

  // Extract data for 2D extraction (bilinear interpolation)
  template <typename T, typename R>
  void operator()(const T &obDat1, const R &obDat2) {
    interpolator.extract(getElement(obDat1), getElement(obDat2));
  }

  // Extract data for 3D extraction (trilinear interpolation)
  template <typename T, typename R, typename U>
  void operator()(const T &obDat1, const R &obDat2, const U &obDat3) {
    interpolator.extract(getFloatElement(obDat1), getFloatElement(obDat2), getFloatElement(obDat3));
  }

  DataExtractor<ExtractedValue> &interpolator;
  size_t iloc;
  int ivar;
};


// -----------------------------------------------------------------------------
template <typename T>
void DrawValueFromFile<T>::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<T> & out) const {
  DataExtractor<T> interpolator{fpath_, options_.group};

  // Channel number handling - do we only process the channels specified in the
  // options, or leave it up to useChannelList variables?
  bool extract_channels = channels_.size() > 0;
  for (auto const & useChannel : useChannelList_) {
    if (useChannel.second) {
      extract_channels = false;
      break;
    }
  }
  if (extract_channels)
    interpolator.scheduleSort("MetaData/sensorChannelNumber", InterpMethod::EXACT,
                              ExtrapolationMode::ERROR, EquidistantChoice::FIRST,
                              CoordinateTransformation::NONE);

  ObData obData;
  for (size_t ind=0; ind < allvars_.size(); ind++) {
    oops::Log::debug() << "Extracting " << allvars_[ind].variable() <<
      " from the obsSpace" << std::endl;

    const std::string varName = allvars_[ind].fullName();
    interpolator.scheduleSort(varName, interpMethod_.at(varName), extrapMode_.at(varName),
                              equidistantChoice_.at(varName),
                              coordinateTransformation_.at(varName));
    if (useChannelList_.at(varName)) {
      switch (in.dtype(allvars_[ind])) {
        case ioda::ObsDtype::Integer:
          updateObData<int>(in, allvars_[ind], obData, channels_);
          break;
        case ioda::ObsDtype::String:
          updateObData<std::string>(in, allvars_[ind], obData, channels_);
          break;
        case ioda::ObsDtype::Float:
          updateObData<float>(in, allvars_[ind], obData, channels_);
          break;
        case ioda::ObsDtype::DateTime:
          updateObDataDateTime(in, allvars_[ind], obData);
          break;
        default:
          throw eckit::UserError("Data type not yet handled.", Here());
      }
    } else {
      switch (in.dtype(allvars_[ind])) {
        case ioda::ObsDtype::Integer:
          updateObData<int>(in, allvars_[ind], obData);
          break;
        case ioda::ObsDtype::String:
          updateObData<std::string>(in, allvars_[ind], obData);
          break;
        case ioda::ObsDtype::Float:
          updateObData<float>(in, allvars_[ind], obData);
          break;
        case ioda::ObsDtype::DateTime:
          updateObDataDateTime(in, allvars_[ind], obData);
          break;
        default:
          throw eckit::UserError("Data type not yet handled.", Here());
      }
    }
  }
  // Finalise (apply) sort by calling with no arguments.
  interpolator.sort();

  for (size_t jvar = 0; jvar < out.nvars(); ++jvar) {
    for (size_t iloc = 0; iloc < in.nlocs(); ++iloc) {
      try {
        if (extract_channels)
          interpolator.extract(channels_[jvar]);

        // Perform any extraction methods.
        ExtractVisitor<T> visitor(interpolator, iloc);
        for (size_t ind=0; ind < obData.size(); ind++) {
          // 'interpolationMethod' is a copy to avoid a MetOffice CRAY icpc compile failure.
          // See https://github.com/JCSDA-internal/ufo/pull/1419
          ufo::InterpMethod interpolationMethod = interpMethod_.at(obData[ind].first);
          visitor.ivar = jvar;
          if ((interpolationMethod == InterpMethod::TRILINEAR) && (ind == (obData.size()-3))) {
            boost::apply_visitor(visitor,
                                 obData[ind].second, obData[ind+1].second, obData[ind+2].second);
            break;
          } else if ((interpolationMethod == InterpMethod::BILINEAR) &&
                     (ind == (obData.size()-2))) {
            boost::apply_visitor(visitor, obData[ind].second, obData[ind+1].second);
            break;
          } else {
            boost::apply_visitor(visitor, obData[ind].second);
          }
        }
        out[jvar][iloc] = interpolator.getResult();
      } catch (const std::exception &ex) {
        // Print extra information that should help the user debug the problem.
        oops::Log::error() << "ERROR: Value extraction failed.\n";
        oops::Log::error() << "  ObsSpace location: " << iloc << "\n";
        oops::Log::error() << "  Interpolation variables:\n";
        // Print values of the interpolation variables at this location
        PrintVisitor visitor(oops::Log::error(), iloc);
        for (size_t ind = 0; ind < obData.size(); ++ind) {
          // Variable name
          oops::Log::error() << "    - " << obData[ind].first << ": ";
          // Variable value
          boost::apply_visitor(visitor, obData[ind].second);
          oops::Log::error() << '\n';
        }
        oops::Log::error() << "  Error message: " << ex.what() << std::endl;
        // Rethrow the original exception
        throw;
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <typename T>
const ufo::Variables & DrawValueFromFile<T>::requiredVariables() const {
  return allvars_;
}

// -----------------------------------------------------------------------------
// Explicit instantiations
template class DrawValueFromFile<float>;
template class DrawValueFromFile<int>;
template class DrawValueFromFile<std::string>;

}  // namespace ufo
