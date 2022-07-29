#ifndef FIXEDPOINT_OPERATORS_HPP
#define FIXEDPOINT_OPERATORS_HPP

#include "hint.hpp"
#include "operations.hpp"
#include "expression_types.hpp"

namespace archgenlib {

#define BINARY_OPERATOR(NAME, BACKEND_NAME, FRONTEND_NAME)                     \
  using NAME##Op = detail::OperationType<OperationKind::NAME>;                 \
  template <ExpressionType T1, ExpressionType T2>                              \
  using NAME##Expr = BinaryOp<T1, T2, NAME##Op>;                               \
  template <archgenlib::ExprHoldType T1, archgenlib::ExprHoldType T2>          \
  auto FRONTEND_NAME(T1 const &op1, T2 const &op2) {                           \
    return archgenlib::detail::binary_op<archgenlib::NAME##Expr>(op1, op2);    \
  }

#define UNARY_OPERATOR(NAME, BACKEND_NAME, FRONTEND_NAME)                      \
  using NAME##Op = detail::OperationType<OperationKind::NAME>;                 \
  template <ExpressionType T> using NAME##Expr = UnaryOp<T, NAME##Op>;         \
  template <archgenlib::ExprHoldType T> auto FRONTEND_NAME(T const &op1) {     \
    return archgenlib::detail::unary_op<archgenlib::NAME##Expr>(op1);          \
  }

#define CONSTANT_OPERATOR(NAME, BACKEND_NAME, FRONTEND_NAME)                   \
  using NAME##Op = detail::OperationType<OperationKind::NAME>;                 \
  using NAME##Expr = NullaryOp<NAME##Op>;                                      \
  NAME##Expr FRONTEND_NAME{};

#include "operators.def"

} // namespace archgenlib

#endif
