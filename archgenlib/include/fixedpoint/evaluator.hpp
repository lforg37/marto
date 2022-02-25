#ifndef FIXEDPOINT_EVALUATOR_HPP
#define FIXEDPOINT_EVALUATOR_HPP

#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>

#include "expression_types.hpp"

namespace archgenlib {
namespace detail {
class SpecialisationFormatter {
public:
  std::ofstream output;
  SpecialisationFormatter();
  ~SpecialisationFormatter();
};
SpecialisationFormatter &getFormatter();
} // namespace detail

template <ExpressionType ET, std::int32_t prec> class Evaluator {
public:
  Evaluator() {
    auto &formatter = detail::getFormatter();
    using introspector_t = detail::Introspecter<ET>;
    if (formatter.output) {
      auto expr_name = introspector_t::get_expression_full_name();
      formatter.output << "template<>\n"
                       << "struct archgenlib::Evaluator<" << expr_name << ", "
                       << prec << "> {\n"
                       << "  auto evaluate(" << expr_name
                       << " const & expr) {\n"
                       << "    return 42;\n"
                       << "  }\n"
                       << "};\n";
    }
  }
  auto evaluate(ET const &expr) { return 57; }
};
} // namespace archgenlib

#endif
