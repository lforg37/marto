#ifndef FIXEDPOINT_SUM_HPP
#define FIXEDPOINT_SUM_HPP

#include "expression_types.hpp"
#include "operations.hpp"
#include <string_view>

namespace archgenlib {

using SumOp = detail::OperationType<OperationKind::ADD>;

template <ExpressionType T1, ExpressionType T2>
using SumExpr = BinaryOp<T1, T2, SumOp>;

using SubOp = detail::OperationType<OperationKind::SUB>;

template <ExpressionType T1, ExpressionType T2>
using SubExpr = BinaryOp<T1, T2, SubOp>;
} // namespace archgenlib

#endif
