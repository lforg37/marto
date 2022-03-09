#ifndef FIXEDPOINT_EVALUATOR_HPP
#define FIXEDPOINT_EVALUATOR_HPP

#include <cstdint>
#include <fstream>
#include <sollya.h>
#include <sstream>
#include <string>

#include "expression_types.hpp"
#include "operators/value_getter.hpp"
#include "operators/tabled_expr.hpp"
#include "runtime/expression_tree.hpp"
#include "runtime/output_formatter.hpp"
#include "runtime/sollya_handler.hpp"
#include "runtime/sollya_operation.hpp"
#include "runtime/sollya_fix_function.hpp"

namespace archgenlib {

template <ExpressionType ET, std::int32_t prec> class Evaluator {
public:
  Evaluator() {
    auto &formatter = detail::getFormatter();
    ExpressionRTRepr erepr{ExprTypeHolder<ET>{}};
    auto getter = ExprLeafGetter::getOperator(
        erepr.symbol_table.at(0).path_from_root);
    auto l = erepr.get_singlevar_dominants();
    SollyaHandler topnode{sollyaFunctionFromNode(*erepr.root)};
    sollya_lib_printf("sollya_repr: %b\n", static_cast<sollya_obj_t>(topnode));
    FPDimRTRepr repr{5, prec, true};
    SollyaFunction sf{topnode, repr};
    auto reprvec = sf.faithful_at_weight(prec);
    if (formatter.output) {
      auto expr_name = detail::type_name<ET>();
      formatter.output << "template<>\n"
                       << "struct archgenlib::Evaluator<" << expr_name << ", "
                       << prec << "> {\n"
                       << "  auto evaluate(" << expr_name
                       << " const & expr) {\n"
                       << "    " << emit_tabled_expr("table", reprvec)
                       << "    return table(static_cast<int>(" << getter("expr")
                       << ".value().unravel()));\n"
                       << "  }\n"
                       << "};\n";
    }
  }
  auto evaluate(ET const &expr) { return 57; }
};
} // namespace archgenlib

#endif
