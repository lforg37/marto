#ifndef FIXEDPOINT_SUM_HPP
#define FIXEDPOINT_SUM_HPP

#include "expression_types.hpp"
#include "operations.hpp"

namespace archgenlib {

using SumOp = detail::OperationType<OperationKind::ADD>;

template <ExpressionType T1, ExpressionType T2>
using SumExpr = BinaryOp<T1, T2, SumOp>;

using SubOp = detail::OperationType<OperationKind::SUB>;

template <ExpressionType T1, ExpressionType T2>
using SubExpr = BinaryOp<T1, T2, SubOp>;

using NegOp = detail::OperationType<OperationKind::NEG>;

template<ExpressionType ET>
using NegExpr = UnaryOp<ET, NegOp>;
} // namespace archgenlib

#endif
