/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/utils/parameters/ParameterBase.h"
#include "ufo/utils/parameters/Parameters.h"

namespace ufo {

ParameterBase::ParameterBase(Parameters *parent) {
  if (parent) {
    parent->registerChild(*this);
  }
}

}  // namespace ufo
