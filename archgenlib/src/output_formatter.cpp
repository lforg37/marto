#include "fixedpoint/evaluator.hpp"
#include <sstream>

namespace archgenlib {
namespace detail {

std::string doReturn(std::string_view expr) {
  std::stringstream ss("return ");
  ss << expr << ";";
  return ss.str();
}

SpecialisationFormatter &getFormatter() {
  static SpecialisationFormatter formatter{};
  return formatter;
}

SpecialisationFormatter::SpecialisationFormatter()
    : output{"specialization_header.hpp"} {
  if (output) {
    output << "#ifndef SPECIALIZATION_HEADER_HPP\n#define "
              "SPECIALIZATION_HEADER_HPP\n";
  }
}

SpecialisationFormatter::~SpecialisationFormatter() {
  if (output) {
    output << "#endif\n";
  }
}

} // namespace detail
} // namespace archgenlib
