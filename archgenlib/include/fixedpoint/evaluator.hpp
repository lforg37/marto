#ifndef FIXEDPOINT_EVALUATOR_HPP
#define FIXEDPOINT_EVALUATOR_HPP

#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>

#include "expression_rt_tree.hpp"
#include "expression_types.hpp"
#include "operators/value_getter.hpp"
#include "output_formatter.hpp"
namespace archgenlib {

template <ExpressionType ET, std::int32_t prec> class Evaluator {
public:
  Evaluator() {
    auto &formatter = detail::getFormatter();
    ExpressionRTRepr erepr{};
    erepr.template do_init<ET>();
    auto getter = ExprLeafGetter::getOperator(
        erepr.symbol_table.value()[0].path_from_root);
    auto l = erepr.get_singlevar_dominants();
    if (formatter.output) {
      auto expr_name = detail::type_name<ET>();
      formatter.output << "template<>\n"
                       << "struct archgenlib::Evaluator<" << expr_name << ", "
                       << prec << "> {\n"
                       << "  auto evaluate(" << expr_name
                       << " const & expr) {\n"
                       << "    return static_cast<int>(" << getter("expr")
                       << ".value().unravel()) + " << l.size() << ";\n"
                       << "  }\n"
                       << "};\n";
    }
  }
  auto evaluate(ET const &expr) { return 57; }
};
} // namespace archgenlib

#endif
