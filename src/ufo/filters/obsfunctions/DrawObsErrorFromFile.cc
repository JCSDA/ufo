/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/exception/Exceptions.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/ObsSpaces.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "ufo/filters/obsfunctions/DrawObsErrorFromFile.h"
#include "ufo/ObsTraits.h"


namespace ufo {

constexpr char InterpMethodParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<InterpMethod>
  InterpMethodParameterTraitsHelper::namedValues[];

static ObsFunctionMaker<DrawObsErrorFromFile>
         makerNetCDF_("DrawObsErrorFromFile");


// -----------------------------------------------------------------------------
DrawObsErrorFromFile::DrawObsErrorFromFile(const eckit::LocalConfiguration config)
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

    fpath_ = options_.fpath.value();
    allvars_  = Variables(interpSubConfs);
}


// -----------------------------------------------------------------------------
DrawObsErrorFromFile::~DrawObsErrorFromFile() {}


class ExtractVisitor : public boost::static_visitor<void> {
 public:
  ExtractVisitor(ObsErrorNetCDFHandler &ncHandler,
                 const size_t &iloc) :
    ncHandler(ncHandler), iloc(iloc) {}

  template <typename T>
  void operator()(const std::vector<T> &obDat) {
    auto obVal = obDat[iloc];
    ncHandler.extract(obVal);
  }

  ObsErrorNetCDFHandler &ncHandler;
  const size_t &iloc;
};


// -----------------------------------------------------------------------------
void DrawObsErrorFromFile::compute(const ObsFilterData & in,
                                      ioda::ObsDataVector<float> & out) const {
  const float missing = util::missingValue(missing);
  ASSERT(out.nvars() == 1);

  ObsErrorNetCDFHandler ncHandler{fpath_};

  ObData obData;
  for (size_t ind=0; ind < allvars_.size(); ind++) {
    oops::Log::debug() << "Extracting " << allvars_[ind].variable() <<
      " from the obsSpace" << std::endl;

    const std::string &varName = get_full_name(allvars_[ind]);
    const InterpMethod &interpolationMethod = interpMethod_.at(varName);
    ncHandler.scheduleSort(varName, interpolationMethod);
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
  ncHandler.sort();

  int jvar = 0;
  if (out.nvars() != 1) {
    throw eckit::NotImplemented("Multichannel variables not yet supported.", Here());
  }
  for (size_t iloc = 0; iloc < in.nlocs(); ++iloc) {
    // Perform any extraction methods (exact, nearest and linear interp.)
    for (auto &od : obData) {
      ExtractVisitor visitor(ncHandler, iloc);
      boost::apply_visitor(visitor, od.second);
    }
    out[static_cast<size_t>(jvar)][iloc] = ncHandler.getObsErrorValue();
  }
}

// -----------------------------------------------------------------------------
const ufo::Variables & DrawObsErrorFromFile::requiredVariables() const {
  return allvars_;
}


// -----------------------------------------------------------------------------

}  // namespace ufo
