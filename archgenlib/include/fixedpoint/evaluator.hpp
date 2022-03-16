#ifndef FIXEDPOINT_EVALUATOR_HPP
#define FIXEDPOINT_EVALUATOR_HPP

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sollya.h>
#include <sstream>
#include <string>

#include "expression_types.hpp"
#include "runtime/expression_tree.hpp"
#include "runtime/operator_builder/affectation.hpp"
#include "runtime/operator_builder/table_builder.hpp"
#include "runtime/operator_builder/value_getter.hpp"
#include "runtime/output_formatter.hpp"
#include "runtime/sollya_fix_function.hpp"
#include "runtime/sollya_handler.hpp"
#include "runtime/sollya_operation.hpp"

namespace archgenlib {

template <ExpressionType ET, std::int32_t prec> class Evaluator {
public:
  Evaluator() {
    auto &formatter = detail::getFormatter();
    ExpressionRTRepr erepr{ExprTypeHolder<ET>{}};
    auto getter =
        ExprLeafGetter{erepr.symbol_table.at(0).path_from_root, "expr"};
    auto freevar = Affectation{"input0", getter};
    auto l = erepr.get_singlevar_dominants();
    auto it = find(l.begin(), l.end(), &erepr.root.value());
    if (it == l.end()) {
      std::cerr << "As of now, only function of one free variable are handled"
                << std::endl;
      std::abort();
    }
    SollyaHandler topnode{sollyaFunctionFromNode(*erepr.root)};
    sollya_lib_printf("sollya_repr: %b\n", static_cast<sollya_obj_t>(topnode));
    FPDimRTRepr repr = erepr.symbol_table.at(0).description.dim;
    SollyaFunction sf{topnode, repr};
    auto reprvec = sf.faithful_at_weight(prec);
    TableBuilder optable{reprvec, repr, {sf.msb_output, prec, sf.signed_output}};
    Affectation tableLambda{"table_func", optable};
    if (formatter.output) {
      auto expr_name = detail::type_name<ET>();
      formatter.output << "template<>\n"
                       << "struct archgenlib::Evaluator<" << expr_name << ", "
                       << prec << "> {\n"
                       << "  auto evaluate(" << expr_name
                       << " const & expr) {\n"
                       << "    " << freevar
                       << "    " << tableLambda
                       << "    " << "return " << tableLambda.get_name() << "("
                       << "static_cast<std::size_t>(" << freevar.get_name() 
                       << ".value()));\n"
                       << "  }\n"
                       << "};\n";
    }
  }
  auto evaluate(ET const &expr) { return FPNumber<FPDim<4, -7, false>>{42}; }
};
} // namespace archgenlib

#endif
