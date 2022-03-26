#ifndef FIXEDPOINT_CONSTANTS_HPP
#define FIXEDPOINT_CONSTANTS_HPP

#include "expression_types.hpp"
#include "operations.hpp"

namespace archgenlib {

using PiOp = detail::OperationType<OperationKind::PI>;

using PiExpr = NullaryOp<PiOp>;

PiExpr pi{};

}

#endif
