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
#include "oops/util/Logger.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/DrawValueFromFile.h"

namespace ufo {

namespace {

// -----------------------------------------------------------------------------
typedef std::vector<std::pair<std::string, boost::variant<std::vector<int>,
                                                          std::vector<float>,
                                                          std::vector<std::string>>
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
    for (auto intParam = interpolationParameters.begin();
         intParam != interpolationParameters.end(); ++intParam) {
      const ufo::InterpMethod & method = intParam->method.value();
      if (method == InterpMethod::BILINEAR) {
        nlin++;
      } else if (nlin > 0) {
        throw eckit::UserError("Bilinear interpolation can only be supplied as the final two "
                               "arguments.", Here());
      } else if ((method == InterpMethod::LINEAR) &&
                 (intParam + 1 != interpolationParameters.end())) {
        throw eckit::UserError("Linear interpolation can only be supplied as the very last "
                               "argument.", Here());
      }
      const std::string varName = ioda::convertV1PathToV2Path(intParam->name.value());
      interpSubConfs.push_back(intParam->toConfiguration());
      interpMethod_[varName] = method;
      extrapMode_[varName] = intParam->extrapMode.value();
      equidistantChoice_[varName] = intParam->equidistanceChoice.value();
    }
    if (nlin > 0 && nlin != 2) {
      throw eckit::UserError("Bilinear interpolation requires two variables.", Here());
    }
    // Get channels from options
    if (options_.chlist.value() != boost::none) {
        std::set<int> channels = options_.chlist.value().get();
        channels_ = {std::make_move_iterator(std::begin(channels)),
                     std::make_move_iterator(std::end(channels))};
    }
    fpath_ = options_.fpath.value();
    allvars_  = Variables(interpSubConfs);
}


// -----------------------------------------------------------------------------
template <typename ExtractedValue>
class ExtractVisitor : public boost::static_visitor<void> {
 public:
  ExtractVisitor(DataExtractor<ExtractedValue> &interpolator, size_t iloc) :
    interpolator(interpolator), iloc(iloc) {}

  template <typename T>
  void operator()(const std::vector<T> &obDat) {
    interpolator.extract(obDat[iloc]);
  }

  template <typename T, typename R>
  void operator()(const std::vector<T> &obDat1, const std::vector<R> &obDat2) {
    interpolator.extract(obDat1[iloc], obDat2[iloc]);
  }

  DataExtractor<ExtractedValue> &interpolator;
  size_t iloc;
};


// -----------------------------------------------------------------------------
template <typename T>
void DrawValueFromFile<T>::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<T> & out) const {
  DataExtractor<T> interpolator{fpath_, options_.group};

  // Channel number handling
  if (options_.chlist.value() != boost::none)
    interpolator.scheduleSort("channel_number@MetaData", InterpMethod::EXACT,
                              ExtrapolationMode::ERROR, EquidistantChoice::FIRST);

  ObData obData;
  for (size_t ind=0; ind < allvars_.size(); ind++) {
    oops::Log::debug() << "Extracting " << allvars_[ind].variable() <<
      " from the obsSpace" << std::endl;

    const std::string varName = allvars_[ind].fullName();
    interpolator.scheduleSort(varName, interpMethod_.at(varName), extrapMode_.at(varName),
                              equidistantChoice_.at(varName));
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
  // Finalise (apply) sort by calling with no arguments.
  interpolator.sort();

  for (size_t jvar = 0; jvar < out.nvars(); ++jvar) {
    for (size_t iloc = 0; iloc < in.nlocs(); ++iloc) {
      try {
        if (options_.chlist.value() != boost::none)
          interpolator.extract(channels_[jvar]);

        // Perform any extraction methods.
        ExtractVisitor<T> visitor(interpolator, iloc);
        for (size_t ind=0; ind < obData.size(); ind++) {
          // 'interpolationMethod' is a copy to avoid a MetOffice CRAY icpc compile failure.
          // See https://github.com/JCSDA-internal/ufo/pull/1419
          ufo::InterpMethod interpolationMethod = interpMethod_.at(obData[ind].first);
          if ((interpolationMethod == InterpMethod::BILINEAR) && (ind == (obData.size()-2))) {
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
