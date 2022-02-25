#ifndef FIXEDPOINT_EXPRESSION_TYPES_HPP
#define FIXEDPOINT_EXPRESSION_TYPES_HPP

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string_view>
#include <type_traits>
#include <utility>

namespace archgenlib {

namespace detail {
template <typename T> constexpr bool is_integral_cst = false;

template <std::integral T, T value>
constexpr bool is_integral_cst<std::integral_constant<T, value>> = true;

template <typename T>
concept IntegralConstantType = is_integral_cst<T>;
} // namespace detail

template <typename T>
concept ExpressionType = requires {
  {T::constant};
};

template <ExpressionType op1type, ExpressionType op2type, typename ChildClass>
struct BinaryOp {
public:
  using arg_storage = std::pair<op1type const &, op2type const &>;
  static constexpr bool constant = op1type::constant || op2type::constant;
  BinaryOp(op1type const &op1, op2type const &op2) : args{op1, op2} {}

private:
  arg_storage args;
};

template <ExpressionType T, typename ChildClass> class UnaryOp {
public:
  static constexpr bool constant = T::constant;

private:
  T &op;
};

template <std::integral T> class Variable {
private:
  T const & value;

public:
  Variable(T const & val) : value{val} {}
  static constexpr bool constant = false;
};

template <detail::IntegralConstantType Integral> class Constant {
private:
  using type = Integral;

public:
  static constexpr bool constant = true;
};

namespace detail {
template <ExpressionType> class Introspecter;

template <ExpressionType LeftOpT, ExpressionType RightOpT, typename ChildClass>
class Introspecter<BinaryOp<LeftOpT, RightOpT, ChildClass>> {
  using left_op_introspector = Introspecter<LeftOpT>;
  using right_op_introspector = Introspecter<RightOpT>;

public:
  static auto get_expression_full_name() {
    std::stringstream ss;
    ss << "archgenlib::BinaryOp<"
       << left_op_introspector::get_expression_full_name() << ", "
       << right_op_introspector::get_expression_full_name() << ", "
       << ChildClass::FQN << ">";
    return ss.str();
  }
};

template <ExpressionType OpT, typename ChildClass>
class Introspecter<UnaryOp<OpT, ChildClass>> {
  using op_introspector_t = Introspecter<OpT>;

public:
  static auto get_expression_full_name() {
    std::stringstream ss;
    ss << ChildClass::FQN << "<"
       << op_introspector_t::get_expression_full_name() << ">";
    return ss.str();
  }
};

template <std::integral T> class Introspecter<Variable<T>> {
public:
  static auto get_expression_full_name() {
    // TODO: replace int by real representation
    constexpr std::string_view ret{"archgenlib::Variable<int>"};
    return ret;
  }
};

template <std::integral T, T value>
class Introspecter<Constant<std::integral_constant<T, value>>> {
public:
  static auto get_expression_full_name() {
    // TODO: replace int by real representation
    std::stringstream ss;
    ss << "archgenlib::Constant<std::integral_constant<int, " << value << ">>";
    return ss.str();
  }
};

} // namespace detail

} // namespace archgenlib

#endif // EXPRESSION_HPP
