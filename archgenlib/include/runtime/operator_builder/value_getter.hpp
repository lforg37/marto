#ifndef RUNTIME_OPERATOR_BUILDER_VALUE_GETTER_HPP
#define RUNTIME_OPERATOR_BUILDER_VALUE_GETTER_HPP

#include "runtime/expression_tree.hpp"
#include <ostream>

namespace archgenlib {

namespace detail {
std::ostream& print_path(std::ostream& os, ExpressionRTRepr::path_t const & path);
}

template<typename SubExprType>
struct ExprLeafGetter {
  ExpressionRTRepr::path_t const & path;
  SubExprType const & subexpr;
  friend std::ostream& operator<<(std::ostream& os, ExprLeafGetter const & elg) {
    os << elg.subexpr;
    return detail::print_path(os, elg.path);
  }
  ExprLeafGetter(ExpressionRTRepr::path_t const & path, SubExprType const & subexpr):path{path}, subexpr{subexpr}
  {}
};
} // namespace archgenlib

#endif
