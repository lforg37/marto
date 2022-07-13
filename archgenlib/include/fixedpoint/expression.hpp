#ifndef FIXEDPOINT_EXPRESSION_HPP
#define FIXEDPOINT_EXPRESSION_HPP

#include "fixedpoint/expression_types.hpp"
#include "functions.hpp"
#include "sum.hpp"
#include "multiplication.hpp"


#include "evaluator.hpp"

template <archgenlib::ExprHoldType T1, archgenlib::ExprHoldType T2>
auto operator+(T1 const &op1, T2 const &op2) {
  return archgenlib::detail::binary_op<archgenlib::SumExpr>(op1, op2);
}

template <archgenlib::ExprHoldType T1, archgenlib::ExprHoldType T2>
auto operator-(T1 const &op1, T2 const &op2) {
  return archgenlib::detail::binary_op<archgenlib::SubExpr>(op1, op2);
}

template <archgenlib::ExpressionType ET> auto operator-(ET const &operand) {
  return archgenlib::detail::unary_op<archgenlib::NegExpr>(operand);
}

template <archgenlib::ExprHoldType T1, archgenlib::ExprHoldType T2>
auto operator*(T1 const &op1, T2 const &op2) {
  return archgenlib::detail::binary_op<archgenlib::MulExpr>(op1, op2);
}


template <archgenlib::ExprHoldType T1, archgenlib::ExprHoldType T2>
auto operator/(T1 const &op1, T2 const &op2) {
  return archgenlib::detail::binary_op<archgenlib::DivExpr>(op1, op2);
}

namespace archgenlib {
template <ExpressionType ET, std::int32_t prec>
static Evaluator<ET, prec> _evaluator{};

template <std::int32_t prec, ExprHoldType ET> auto evaluate(ET const &val) {
  auto holder = detail::get_holder(val);
  using h_t = decltype(holder);
  return _evaluator<typename h_t::expression_t, prec>.evaluate(holder.expression);
}
} // namespace archgenlib

#endif
