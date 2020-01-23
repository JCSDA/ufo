/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_ARROWPROXY_H_
#define UFO_UTILS_ARROWPROXY_H_

namespace ufo {

/// \brief Utility class used in overloads of operator-> in forward iterators.
///
/// See e.g. A. O'Dwyer's blog, https://quuxplusone.github.io/blog/2019/02/06/arrow-proxy.
template<class T>
class ArrowProxy {
 public:
  explicit ArrowProxy(const T &ref) :
    ref_(ref)
  {}

  T *operator->() {
    return &ref_;
  }

 private:
  T ref_;
};

}  // namespace ufo

#endif  // UFO_UTILS_ARROWPROXY_H_
