#ifndef FIXED_POINT_OPERATIONS_HPP
#define FIXED_POINT_OPERATIONS_HPP

namespace archgenlib {
enum struct OperationKind { ADD, SUB, MUL, DIV, SIN, LOG2, LOG, ABS, NEG, POW, PI };

namespace detail {
constexpr bool isBinaryOpKind(OperationKind OK) {
  switch (OK) {
  case OperationKind::ADD:
  case OperationKind::SUB:
  case OperationKind::MUL:
  case OperationKind::DIV:
  case OperationKind::POW:
    return true;
  case OperationKind::NEG:
  case OperationKind::SIN:
  case OperationKind::LOG2:
  case OperationKind::LOG:
  case OperationKind::ABS:
  case OperationKind::PI:
    return false;
  default:
    __builtin_unreachable();
  }
}

constexpr bool isNullaryOpKind(OperationKind OK) {
  switch (OK) {
  case OperationKind::ADD:
  case OperationKind::SUB:
  case OperationKind::MUL:
  case OperationKind::DIV:
  case OperationKind::LOG:
  case OperationKind::POW:
  case OperationKind::NEG:
  case OperationKind::SIN:
  case OperationKind::LOG2:
  case OperationKind::ABS:
    return false;
  case OperationKind::PI:
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
