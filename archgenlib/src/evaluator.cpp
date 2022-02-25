#include "fixedpoint/evaluator.hpp"

namespace archgenlib {
namespace detail {
SpecialisationFormatter &getFormatter() {
  static SpecialisationFormatter formatter{};
  return formatter;
}

SpecialisationFormatter::SpecialisationFormatter()
    : output{"specialization_header.hpp"} {
  if (output) {
    output << "#ifndef SPECIALIzATION_HEADER_HPP\n#define "
              "SPECIALIzATION_HEADER_HPP\n";
  }
}

SpecialisationFormatter::~SpecialisationFormatter() {
  if (output) {
    output << "#endif\n";
  }
}

} // namespace detail
} // namespace archgenlib