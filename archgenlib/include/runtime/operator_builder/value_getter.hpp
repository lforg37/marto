#ifndef RUNTIME_OPERATOR_BUILDER_VALUE_GETTER_HPP
#define RUNTIME_OPERATOR_BUILDER_VALUE_GETTER_HPP

#include "runtime/expression_tree.hpp"
#include "runtime/operator_builder/operator.hpp"

namespace archgenlib {

detail::cpp_expr_ptr get_path(ExpressionRTRepr::path_t const &path,
                              detail::cpp_expr_ptr);
} // namespace archgenlib

#endif
