#ifndef FIXEDPOINT_EVALUATOR_HPP
#define FIXEDPOINT_EVALUATOR_HPP

#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>

#include "expression_types.hpp"
#include "expression_rt_tree.hpp"
#include "output_formatter.hpp"
namespace archgenlib {

template <ExpressionType ET, std::int32_t prec> class Evaluator {
public:
  Evaluator() {
    auto &formatter = detail::getFormatter();
    ExpressionRTRepr erepr{};
    erepr.template do_init<ET>();
    if (formatter.output) {
      auto expr_name = detail::type_name<ET>();
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
