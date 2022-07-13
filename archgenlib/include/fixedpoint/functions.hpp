#ifndef FIXEDPOINT_FUNCTIONS_HPP
#define FIXEDPOINT_FUNCTIONS_HPP

#include "expression_types.hpp"
#include "operations.hpp"

namespace archgenlib {
using SinOp = detail::OperationType<OperationKind::SIN>;

template<ExpressionType ET>
using SinExpr = UnaryOp<ET, SinOp>;

template<ExprHoldType ET>
auto sin(ET const & operand) {
  return detail::unary_op<SinExpr>(operand);
}

using Log2Op = detail::OperationType<OperationKind::LOG2>;

template<ExpressionType ET>
using Log2Expr = UnaryOp<ET, Log2Op>;

template<ExprHoldType ET>
auto log2(ET const & operand) {
  return detail::unary_op<Log2Expr>(operand);
}

using LogOp = detail::OperationType<OperationKind::LOG>;

template<ExpressionType ET>
using LogExpr = UnaryOp<ET, LogOp>;

template<ExprHoldType ET>
auto log(ET const & operand) {
  return detail::unary_op<LogExpr>(operand);
}

using AbsOp = detail::OperationType<OperationKind::ABS>;

template<ExpressionType ET>
using AbsExpr = UnaryOp<ET, AbsOp>;

template<ExprHoldType ET>
auto abs(ET const & operand) {
  return detail::unary_op<AbsExpr>(operand);
}

using PowOp = detail::OperationType<OperationKind::POW>;

template<ExpressionType ET, ExpressionType ET2>
using PowExpr = BinaryOp<ET, ET2, PowOp>;

template <ExprHoldType T1, ExprHoldType T2>
auto pow(T1 const &op1, T2 const &op2) {
  return detail::binary_op<PowExpr>(op1, op2);
}

}

#endif
