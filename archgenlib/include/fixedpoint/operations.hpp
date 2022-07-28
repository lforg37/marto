#ifndef FIXED_POINT_OPERATIONS_HPP
#define FIXED_POINT_OPERATIONS_HPP

namespace archgenlib {
enum struct OperationKind {
#define OPERATOR(NAME, BACKEND_NAME, FUNCTION_NAME) NAME,
#include "operators.def"
};

namespace detail {
constexpr bool isBinaryOpKind(OperationKind OK) {
  switch (OK) {
#define BINARY_OPERATOR(NAME, BACKEND_NAME, FUNCTION_NAME)                     \
  case OperationKind::NAME:
#include "operators.def"
    return true;
#define UNARY_OPERATOR(NAME, BACKEND_NAME, FUNCTION_NAME)                      \
  case OperationKind::NAME:
#define CONSTANT_OPERATOR(NAME, BACKEND_NAME, FUNCTION_NAME)                   \
  case OperationKind::NAME:
#include "operators.def"
    return false;
  default:
    __builtin_unreachable();
  }
}

constexpr bool isNullaryOpKind(OperationKind OK) {
  switch (OK) {
#define BINARY_OPERATOR(NAME, BACKEND_NAME, FUNCTION_NAME)                     \
  case OperationKind::NAME:
#define UNARY_OPERATOR(NAME, BACKEND_NAME, FUNCTION_NAME)                      \
  case OperationKind::NAME:
#include "operators.def"
    return false;
#define CONSTANT_OPERATOR(NAME, BACKEND_NAME, FUNCTION_NAME)                   \
  case OperationKind::NAME:
#include "operators.def"
    return true;
  default:
    __builtin_unreachable();
  }
}

constexpr bool isUnaryOpKind(OperationKind OK) {
  return !isBinaryOpKind(OK) && !isNullaryOpKind(OK);
}

template <OperationKind OK> struct OperationType {
  static constexpr bool is_binary = isBinaryOpKind(OK);
  static constexpr bool is_unary = isUnaryOpKind(OK);
  static constexpr bool is_nullary = isNullaryOpKind(OK);
  static constexpr OperationKind operation_kind = OK;
};

template <typename T> constexpr bool _is_binary_op = false;
template <typename T> constexpr bool _is_unary_op = false;
template <typename T> constexpr bool _is_nullary_op = false;

template <OperationKind OK>
constexpr bool
    _is_binary_op<detail::OperationType<OK>> = detail::isBinaryOpKind(OK);

template <OperationKind OK>
constexpr bool
    _is_unary_op<detail::OperationType<OK>> = detail::isUnaryOpKind(OK);    

template <OperationKind OK>
constexpr bool
    _is_nullary_op<detail::OperationType<OK>> = detail::isNullaryOpKind(OK);    

} // namespace detail

template<typename T>
concept BinaryOpType = detail::_is_binary_op<T>;

template<typename T>
concept UnaryOpType = detail::_is_unary_op<T>;

template<typename T>
concept NullaryOpType = detail::_is_nullary_op<T>;

}; // namespace archgenlib

#endif
