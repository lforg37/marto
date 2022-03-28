#ifndef FIXEDPOINT_FUNCTIONS_HPP
#define FIXEDPOINT_FUNCTIONS_HPP

#include "expression_types.hpp"
#include "operations.hpp"

namespace archgenlib {
using SinOp = detail::OperationType<OperationKind::SIN>;

template<ExpressionType ET>
using SinExpr = UnaryOp<ET, SinOp>;

template<ExpressionType ET>
SinExpr<ET> sin(ET const & operand) {
  return {operand};
}

using Log2Op = detail::OperationType<OperationKind::LOG2>;

template<ExpressionType ET>
using Log2Expr = UnaryOp<ET, Log2Op>;

template<ExpressionType ET>
Log2Expr<ET> log2(ET const & operand) {
  return {operand};
}

using LogOp = detail::OperationType<OperationKind::LOG>;

template<ExpressionType ET>
using LogExpr = UnaryOp<ET, LogOp>;

template<ExpressionType ET>
LogExpr<ET> log(ET const & operand) {
  return {operand};
}

using AbsOp = detail::OperationType<OperationKind::ABS>;

template<ExpressionType ET>
using AbsExpr = UnaryOp<ET, AbsOp>;

template<ExpressionType ET>
AbsExpr<ET> abs(ET const & operand) {
  return {operand};
}

using PowOp = detail::OperationType<OperationKind::POW>;

template<ExpressionType ET, ExpressionType ET2>
using PowExpr = BinaryOp<ET, ET2, PowOp>;

template<ExpressionType ET, ExpressionType ET2>
PowExpr<ET, ET2> pow(ET const & operand, ET2 const & operand2) {
  return {operand, operand2};
}

}

#endif
