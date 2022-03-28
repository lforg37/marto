#ifndef FIXEDPOINT_EXPRESSION_HPP
#define FIXEDPOINT_EXPRESSION_HPP

#include "functions.hpp"
#include "sum.hpp"
#include "multiplication.hpp"


#include "evaluator.hpp"

template <archgenlib::ExpressionType T1, archgenlib::ExpressionType T2>
auto operator+(T1 const &op1, T2 const &op2) {
  using ret_type = archgenlib::SumExpr<T1, T2>;
  return ret_type{op1, op2};
}

template <archgenlib::ExpressionType T1, archgenlib::ExpressionType T2>
auto operator-(T1 const &op1, T2 const &op2) {
  using ret_type = archgenlib::SubExpr<T1, T2>;
  return ret_type{op1, op2};
}

template <archgenlib::ExpressionType ET> auto operator-(ET const &operand) {
  return archgenlib::NegExpr<ET>(operand);
}

template <archgenlib::ExpressionType T1, archgenlib::ExpressionType T2>
auto operator*(T1 const &op1, T2 const &op2) {
  using ret_type = archgenlib::MulExpr<T1, T2>;
  return ret_type{op1, op2};
}


template <archgenlib::ExpressionType T1, archgenlib::ExpressionType T2>
auto operator/(T1 const &op1, T2 const &op2) {
  using ret_type = archgenlib::DivExpr<T1, T2>;
  return ret_type{op1, op2};
}

namespace archgenlib {
template <ExpressionType ET, std::int32_t prec>
static Evaluator<ET, prec> _evaluator{};

template <std::int32_t prec, ExpressionType ET> auto evaluate(ET const &val) {
  return _evaluator<ET, prec>.evaluate(val);
}
} // namespace archgenlib

#endif
