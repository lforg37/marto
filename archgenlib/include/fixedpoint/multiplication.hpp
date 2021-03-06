#ifndef FIXEDPOINT_MULTIPLICATION_HPP
#define FIXEDPOINT_MULTIPLICATION_HPP

#include "expression_types.hpp"
#include "operations.hpp"

namespace archgenlib {

using MulOp = detail::OperationType<OperationKind::MUL>;

template <ExpressionType T1, ExpressionType T2>
using MulExpr = BinaryOp<T1, T2, MulOp>;


using DivOp = detail::OperationType<OperationKind::DIV>;

template <ExpressionType T1, ExpressionType T2>
using DivExpr = BinaryOp<T1, T2, DivOp>;

} // namespace archgenlib

#endif
