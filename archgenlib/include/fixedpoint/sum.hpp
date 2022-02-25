#ifndef FIXEDPOINT_SUM_HPP
#define FIXEDPOINT_SUM_HPP

#include "expression_types.hpp"
#include <string_view>

namespace archgenlib {
class SumOp {
public:
  static constexpr std::string_view FQN{"archgenlib::SumOp"};
};

template <ExpressionType T1, ExpressionType T2>
using SumExpr = BinaryOp<T1, T2, SumOp>;

class SunOp {
public:
  static constexpr std::string_view FQN{"archgenlib::SubOp"};
};

template <ExpressionType T1, ExpressionType T2>
using SubExpr = BinaryOp<T1, T2, SumOp>;
} // namespace archgenlib

#endif
