#ifndef FIXEDPOINT_EVALUATE_HPP
#define FIXEDPOINT_EVALUATE_HPP

#include "operators.hpp"
#ifndef ARCHGEN_USE_MLIR_PLUGIN
#include "header_gen_evaluator.hpp"
#else
#include "plugin_gen_evaluator.hpp"
#endif

namespace archgenlib {

#ifndef ARCHGEN_USE_MLIR_PLUGIN
template <ExpressionType ET, std::int32_t prec>
static Evaluator<ET, prec> _evaluator{};
#endif

template <typename OutTy, ExprHoldType ET>
#ifndef ARCHGEN_USE_MLIR_PLUGIN
auto
#else
archgenlib::FixedNumber<OutTy>
#endif
evaluate(ET const &val) {
  auto holder = detail::get_holder(val);
  using h_t = decltype(holder);
#ifndef ARCHGEN_USE_MLIR_PLUGIN
  return _evaluator<typename h_t::expression_t, OutTy::lsb_weight>.evaluate(holder.expression);
#else
  return detail::evaluate<archgenlib::FixedNumber<OutTy>, typename h_t::expression_t>(holder.expression.parameters);
#endif
}
} // namespace archgenlib

#endif
