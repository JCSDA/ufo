/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "eckit/exception/Exceptions.h"
#include "ufo/filters/obsfunctions/DrawValueFromFile.h"


namespace ufo {

namespace {

// -----------------------------------------------------------------------------
typedef std::list<std::pair<std::string, boost::variant<std::vector<int>,
                                                        std::vector<float>,
                                                        std::vector<std::string>>
                           >> ObData;


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
    for (auto intParam = interpolationParameters.begin();
         intParam != interpolationParameters.end(); ++intParam) {
      const ufo::InterpMethod & method = intParam->method.value();
      if ((method == InterpMethod::LINEAR) && (intParam + 1 != interpolationParameters.end())) {
        throw eckit::UserError("Linear interpolation can only be supplied as the very last "
                               "argument.", Here());
      }
      interpSubConfs.push_back(intParam->toConfiguration());
      interpMethod_[intParam->name.value()] = method;
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
  ExtractVisitor(DataExtractor<ExtractedValue> &interpolator,
                 const size_t &iloc) :
    interpolator(interpolator), iloc(iloc) {}

  template <typename T>
  void operator()(const std::vector<T> &obDat) {
    auto obVal = obDat[iloc];
    interpolator.extract(obVal);
  }

  DataExtractor<ExtractedValue> &interpolator;
  const size_t &iloc;
};


// -----------------------------------------------------------------------------
template <typename T>
void DrawValueFromFile<T>::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<T> & out) const {
  DataExtractor<T> interpolator{fpath_, options_.group};

  // Channel number handling
  if (options_.chlist.value() != boost::none)
    interpolator.scheduleSort("channel_number@MetaData", InterpMethod::EXACT);

  ObData obData;
  for (size_t ind=0; ind < allvars_.size(); ind++) {
    oops::Log::debug() << "Extracting " << allvars_[ind].variable() <<
      " from the obsSpace" << std::endl;

    const std::string varName = allvars_[ind].fullName();
    const InterpMethod &interpolationMethod = interpMethod_.at(varName);
    interpolator.scheduleSort(varName, interpolationMethod);
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
      if (options_.chlist.value() != boost::none)
        interpolator.extract(channels_[jvar]);

      // Perform any extraction methods (exact, nearest and linear interp.)
      for (auto &od : obData) {
        ExtractVisitor<T> visitor(interpolator, iloc);
        boost::apply_visitor(visitor, od.second);
      }
      out[jvar][iloc] = interpolator.getResult();
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
