#ifndef FIXED_POINT_OPERATIONS_HPP
#define FIXED_POINT_OPERATIONS_HPP

namespace archgenlib {
enum struct OperationKind { ADD, SUB, SIN, NEG };

namespace detail {
template <OperationKind OK> constexpr bool isBinaryOpKind() {
  switch (OK) {
  case OperationKind::ADD:
  case OperationKind::SUB:
    return true;
  case OperationKind::NEG:
  case OperationKind::SIN:
    return false;
  default:
    __builtin_unreachable();
  }
}

template <OperationKind OK> struct OperationType {
  static constexpr bool is_binary = isBinaryOpKind<OK>();
  static constexpr bool is_unary = !isBinaryOpKind<OK>();
  static constexpr OperationKind operation_kind = OK;
};

template <typename T> constexpr bool _is_binary_op = false;
template <typename T> constexpr bool _is_unary_op = false;

template <OperationKind OK>
constexpr bool
    _is_binary_op<detail::OperationType<OK>> = detail::isBinaryOpKind<OK>();

template <OperationKind OK>
constexpr bool
    _is_unary_op<detail::OperationType<OK>> = !detail::isBinaryOpKind<OK>();    
} // namespace detail

template<typename T>
concept BinaryOpType = detail::_is_binary_op<T>;

template<typename T>
concept UnaryOpType = detail::_is_unary_op<T>;

}; // namespace archgenlib

#endif
