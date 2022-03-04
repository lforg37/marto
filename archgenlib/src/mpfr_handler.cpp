#include "runtime/mpfr_handler.hpp"

namespace archgenlib {

MPFRHandler::MPFRHandler(size_t precision) { mpfr_init2(managed, precision); }
MPFRHandler::~MPFRHandler() {
  mpfr_clear(managed);
}
} // namespace archgenlib
