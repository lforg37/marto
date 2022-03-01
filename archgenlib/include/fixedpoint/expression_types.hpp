#ifndef FIXEDPOINT_EXPRESSION_TYPES_HPP
#define FIXEDPOINT_EXPRESSION_TYPES_HPP

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string_view>
#include <type_traits>
#include <utility>

#include "operations.hpp"

namespace archgenlib {

namespace detail {
template <typename T> constexpr bool is_integral_cst = false;

template <std::integral T, T value>
constexpr bool is_integral_cst<std::integral_constant<T, value>> = true;

template <typename T>
concept IntegralConstantType = is_integral_cst<T>;
} // namespace detail

template <typename T>
concept ExpressionType = requires {
                           { T::constant };
                         };

template <ExpressionType op1type, ExpressionType op2type, BinaryOpType OpT>
struct BinaryOp {
public:
  using LeftType = op1type;
  using RightType = op2type;
  using operation_t = OpT;
  static constexpr bool constant = op1type::constant || op2type::constant;
  BinaryOp(op1type const &op1, op2type const &op2) : left{op1}, right{op2} {}
  LeftType const &left;
  RightType const &right;
};

namespace detail {
template <ExpressionType T> constexpr bool _is_binary_expr = false;

template <ExpressionType ET1, ExpressionType ET2, BinaryOpType OT>
constexpr bool _is_binary_expr<BinaryOp<ET1, ET2, OT>> = true;
} // namespace detail

template <typename ET>
concept BinaryExprType = detail::_is_binary_expr<ET>;

template <ExpressionType T, UnaryOpType Op> class UnaryOp {
public:
  static constexpr bool constant = T::constant;
  using operation_t = Op;
  using ChildType = T;

private:
  T &op;
};

namespace detail {
template <ExpressionType T> constexpr bool _is_unary_expr = false;

template <ExpressionType ET, UnaryOpType OT>
constexpr bool _is_unary_expr<UnaryOp<ET, OT>> = true;
} // namespace detail

template <typename ET>
concept UnaryExprType = detail::_is_unary_expr<ET>;

template <std::integral T> class Variable {
public:
  using type = T;
  Variable(T const &val) : value{val} {}
  static constexpr bool constant = false;
  T const &value;
};

namespace detail {
template <typename T> constexpr bool _is_variable_expr = false;
template <std::integral IT>
constexpr bool _is_variable_expr<Variable<IT>> = true;
} // namespace detail

template <typename T>
concept VariableExprType = detail::_is_variable_expr<T>;

template <detail::IntegralConstantType Integral> class Constant {
public:
  using type = Integral;
  typename Integral::type get() { return type::value; }
  static constexpr bool constant = true;
};

namespace detail {
template <typename T> constexpr bool _is_constant_expr = false;
template <detail::IntegralConstantType IT>
constexpr bool _is_constant_expr<Constant<IT>> = true;
} // namespace detail

template <typename T>
concept ConstantExprType = detail::_is_constant_expr<T>;

} // namespace archgenlib

#endif // EXPRESSION_HPP
