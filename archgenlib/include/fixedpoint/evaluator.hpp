#ifndef FIXEDPOINT_EVALUATOR_HPP
#define FIXEDPOINT_EVALUATOR_HPP

#ifndef INCLUDE_GENERATED_HEADER
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sollya.h>
#include <sstream>
#include <string>

#include "runtime/expression_tree.hpp"
#include "runtime/fixfunction_multipartite.hpp"
#include "runtime/operator_builder/multipartite_operator.hpp"
#include "runtime/operator_builder/operator.hpp"
#include "runtime/operator_builder/table_builder.hpp"
#include "runtime/operator_builder/value_getter.hpp"
#include "runtime/output_formatter.hpp"
#include "runtime/sollya_fix_function.hpp"
#include "runtime/sollya_handler.hpp"
#include "runtime/sollya_operation.hpp"
#endif

#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"

namespace archgenlib {

#ifdef INCLUDE_GENERATED_HEADER
constexpr bool has_specialization_header = true;
#else
constexpr bool has_specialization_header = false;
#endif

template <typename ET, std::int32_t prec> class Evaluator {
public:
  Evaluator() {
#ifndef INCLUDE_GENERATED_HEADER
    auto &formatter = detail::getFormatter();
    ExpressionRTRepr erepr{ExprTypeHolder<ET>{}};
    detail::cpp_expr_ptr input_var{new detail::NamedExpression{"expr"}};
    auto& variable_desc = erepr.symbol_table.begin()->second;
    auto input_val =
        get_path(variable_desc.path_from_root, input_var);
    auto freevar = input_val->assign_to("input0");
    auto l = erepr.get_singlevar_dominants();
    auto it = find(l.begin(), l.end(), &erepr.root.value());
    if (it == l.end()) {
      std::cerr << "As of now, only function of one free variable are handled"
                << std::endl;
      std::abort();
    }
    SollyaHandler topnode{sollyaFunctionFromNode(*erepr.root)};
    sollya_lib_printf("sollya_repr: %b\n", static_cast<sollya_obj_t>(topnode));
    FixedFormatRTRepr repr = variable_desc.description.dim;
    SollyaFunction sf{topnode, repr};
    auto reprvec = sf.faithful_at_weight(prec);
    MultipartiteFunction mpf{sf, prec};
    detail::cpp_expr_ptr table_op;
    std::string_view op_name;
    if (!mpf.best_config.has_value()) {
      auto vals = sf.faithful_at_weight(prec);
      TableBuilder tb{repr, {sf.msb_output, prec, sf.signed_output}, vals};
      table_op.reset(new CPPOperator{tb.build_table()});
      op_name = "tabulated_func";
    } else {
      MultipartiteOperator mop{mpf};
      table_op.reset(new CPPOperator{mop.get_operator()});
      op_name = "bipartite_decomposition";
    }

    auto op = table_op->assign_to(op_name);
    auto ret = (*op)(freevar)->to_return();
    if (formatter.output) {
      auto expr_name = detail::type_name<ET>();
      formatter.output << "template<>\n"
                       << "struct archgenlib::Evaluator<" << expr_name << ", "
                       << prec << "> {\n"
                       << "  auto evaluate(" << expr_name
                       << " const & expr) {\n"
                       << "    " << *freevar
                       << "    " << *op
                       << "    " << *ret
                       << "  }\n"
                       << "};\n";
    }
#endif
  }
  auto evaluate(ET const &expr) { return FixedNumber<FixedFormat<4, -7, unsigned>>{42}; }
};
} // namespace archgenlib

#define STRINGIFY(X) STRINGIFY2(X)
#define STRINGIFY2(X) #X

#ifdef INCLUDE_GENERATED_HEADER
#include STRINGIFY(INCLUDE_GENERATED_HEADER)
#endif

#endif
