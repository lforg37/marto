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

}

#endif