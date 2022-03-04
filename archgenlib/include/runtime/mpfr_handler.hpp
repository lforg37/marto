#ifndef RUNTIME_MPFR_HANDLER_HPP
#define RUNTIME_MPFR_HANDLER_HPP

#include <memory>
#include <type_traits>

#include <mpfr.h>

namespace archgenlib {
  class MPFRHandler {
private:
  mpfr_t managed;

public:
  MPFRHandler(size_t precision);
  MPFRHandler(MPFRHandler const &) = delete;
  MPFRHandler(MPFRHandler &&) = delete;
  ~MPFRHandler();
  operator mpfr_t& () { return managed; };
};
}

#endif
