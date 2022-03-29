#include "fixedpoint/evaluator.hpp"
#include <sstream>
#include <cstdlib>

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
    : output_file{std::getenv("ARCHGENLIB_SPECIALIZATION_HEADER_PATH")}, output{output_file} {
  std::cout << "Generating " << output_file << std::endl;
  if (output) {
    output << "#ifndef SPECIALIZATION_HEADER_HPP\n#define "
           << "SPECIALIZATION_HEADER_HPP\n\n"
           << "#include \"fixedpoint/operators.hpp\"\n\n";
  }
}

SpecialisationFormatter::~SpecialisationFormatter() {
  if (output) {
    output << "#endif\n";
  }
}

} // namespace detail
} // namespace archgenlib
