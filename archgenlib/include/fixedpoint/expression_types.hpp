#ifndef FIXEDPOINT_EXPRESSION_TYPES_HPP
#define FIXEDPOINT_EXPRESSION_TYPES_HPP

#include <bits/iterator_concepts.h>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

#include "fixedpoint.hpp"
#include "operations.hpp"

namespace archgenlib {

namespace detail {
template <typename T> constexpr bool is_integral_cst = false;

template <std::integral T, T value>
constexpr bool is_integral_cst<std::integral_constant<T, value>> = true;

template <typename T>
concept IntegralConstantType = is_integral_cst<T>;

template<typename, typename> struct TupleMerger {};

template<typename... T1, typename... T2> struct TupleMerger<std::tuple<T1...>, std::tuple<T2...>> {
  using result_t = std::tuple<T1..., T2...>;
};

template<typename T1, typename T2>
using concat_tuple_t = typename TupleMerger<T1, T2>::result_t;

/**
 * @brief Traits for signaling expression type
 */
template<typename T>
constexpr bool _is_expression_type = false;

} // namespace detail

/**
 * @brief Concept that is true for types that can appear in an expression tree
 * 
 * @tparam T 
 */
template <typename T>
concept ExpressionType = detail::_is_expression_type<T>;

/**
 * @brief Binary operation node of expression tree
 * 
 * @tparam op1type Left subexpression
 * @tparam op2type Right sub expression
 * @tparam OpT kind of the operation
 */
template <ExpressionType op1type, ExpressionType op2type, BinaryOpType OpT>
struct BinaryOp {
public:
  using LeftType = op1type;
  using RightType = op2type;
  using operation_t = OpT;
  using parameters_t = detail::concat_tuple_t<typename LeftType::parameters_t, typename RightType::parameters_t>;
  static constexpr bool constant = op1type::constant || op2type::constant;
  BinaryOp(op1type const &op1, op2type const &op2) : parameters{std::tuple_cat(op1.parameters, op2.parameters)} {}
  const parameters_t parameters;
};

namespace detail {
template <typename T> constexpr bool _is_binary_expr = false;

template <ExpressionType ET1, ExpressionType ET2, BinaryOpType OT>
constexpr bool _is_binary_expr<BinaryOp<ET1, ET2, OT>> = true;
} // namespace detail

template <typename ET>
concept BinaryExprType = detail::_is_binary_expr<ET>;

namespace detail {
template<BinaryExprType BET>
constexpr bool _is_expression_type<BET> = true;
}

/**
 * @brief Concept for unary expression
 * 
 * @tparam T subexpression to which the operation is applied
 * @tparam Op kind of the operation
 */
template <ExpressionType T, UnaryOpType Op> class UnaryOp {
public:
  static constexpr bool constant = T::constant;
  using operation_t = Op;
  using ChildType = T;
  using parameters_t = typename T::parameters_t;
  UnaryOp(T const & op):parameters{op.parameters}{}
  parameters_t const parameters;
};

namespace detail {
template <typename T> constexpr bool _is_unary_expr = false;

template <ExpressionType ET, UnaryOpType OT>
constexpr bool _is_unary_expr<UnaryOp<ET, OT>> = true;
} // namespace detail

/**
 * @brief Concept for unary expression 
 */
template <typename ET>
concept UnaryExprType = detail::_is_unary_expr<ET>;

namespace detail {
template<UnaryExprType UET>
constexpr bool _is_expression_type<UET> = true;
}

struct ConstantExpr {
  static constexpr bool constant = true;
  using parameters_t = std::tuple<>;
  parameters_t parameters;
};

//TODO: change the name
/**
 * @brief Class for irrational (need of infinite precision) constants
 * 
 * @tparam Op constant kind
 */
template<NullaryOpType Op>
class NullaryOp : public ConstantExpr {
  public:
  using operation_t = Op;
  NullaryOp() {}
};

namespace detail {
template <typename T> constexpr bool _is_nullary_op_expr = false;

template <NullaryOpType OT>
constexpr bool _is_nullary_op_expr<NullaryOp<OT>> = true;
} // namespace detail

/**
 * @brief Concept for identifying irrational constant 
 */
template <typename ET>
concept NullaryOpExprType = detail::_is_nullary_op_expr<ET>;

namespace detail {
template<NullaryOpExprType NOET>
constexpr bool _is_expression_type<NOET> = true;
}

template <typename DisambiguationMaker, typename VET>
struct DisambiguationHolder;

/**
 * @brief Represents a variable node in an expression
 * 
 * @tparam T The format of the variable
 * @tparam VariableID the id of the variable in the current expression tree
 */
template <FixedNumberType T, std::size_t VariableID> class Variable {
public:
  using dimension_t = typename T::format_t;
  using type = T;
  using parameters_t = std::tuple<T const &>;
  static constexpr bool constant = false;
  parameters_t const parameters;
private:
  Variable(T const &val) : parameters{std::forward_as_tuple(val)} {}
  template <typename, typename>
  friend struct DisambiguationHolder;
};

namespace detail {
template <typename T> constexpr bool _is_variable_expr = false;
template <FixedNumberType IT, std::size_t VarID>
constexpr bool _is_variable_expr<Variable<IT, VarID>> = true;
template <FixedNumberType IT, std::size_t VarID>
constexpr bool _is_expression_type<Variable<IT, VarID>> = true;
} // namespace detail

/**
 * @brief Concept for Variable in an expression tree 
 */
template <typename T>
concept VariableExprType = ExpressionType<T> && detail::_is_variable_expr<T>;

/**
 * @brief Represents a fixed constant in the expression tree
 * 
 * @tparam ConstType the constant associated to this expression type
 */
template <FixedConstantType ConstType> class Constant : public ConstantExpr {
public:
  using dimension_t = typename ConstType::dimension_t;
  using type = ConstType;
};

namespace detail {
template <typename T> constexpr bool _is_constant_expr = false;
template <FixedConstantType ConstType>
constexpr bool _is_constant_expr<Constant<ConstType>> = true;
} // namespace detail

template <typename T>
concept ConstantExprType = detail::_is_constant_expr<T>;

namespace detail {
template<ConstantExprType CET>
constexpr bool _is_expression_type<CET> = true;
}

/** Disambiguation of variable types
 *
 *  Expression are identified by the type of the Expression Tree (ET).
 *  The leaves of the ET are constant and variables.
 *  Constant type contains the value, so there is no ambiguities on this.
 *  For variables, the value is only known at runtime, so the value cannot
 *  be part of the type.
 *  
 *  In order to be able to identify uniquely variables of an expression sharing the same type, 
 *  each FreeVariable, each expression tree is associated with a disambiguation type
 */
/// 
//-------------------- Handling disambiguation of variable types -----------//

namespace detail {
struct NoAmbiguities{};
}
template<typename ET>
concept UniqueExprType = ExpressionType<ET> && ET::constant;

template<ExpressionType ET, typename DisambiguationMarker>
struct DisambiguationHolder<ET, DisambiguationMarker> {
  using disambiguation_t = DisambiguationMarker;
  using expression_t = ET;
  expression_t const expression;
};

template <VariableExprType VET, typename DisambiguationMarker>
struct DisambiguationHolder<VET, DisambiguationMarker> {
  using disambiguation_t = DisambiguationMarker;
  using expression_t = VET;
  expression_t const expression;
private:
  DisambiguationHolder(typename expression_t::type const & val):expression{val}{}
  template<typename MarkingType, FixedNumberType... FNT, std::size_t... indexes>
  friend auto get_free_variables(FNT const & ... fnt, std::index_sequence<indexes...>);
  template<typename, FixedNumberType FNT>
  friend auto FreeVariable(FNT const & fnt);
};

namespace detail {
template <typename>
constexpr bool _is_holder = false;

template<ExpressionType ET, typename T>
constexpr bool _is_holder<DisambiguationHolder<ET, T>> = true;
}

template<typename T>
concept DisambiguationHType = detail::_is_holder<T>;

template <UniqueExprType ET>
using NonAmbiguousHolder = DisambiguationHolder<ET, detail::NoAmbiguities>;

template<VariableExprType VET, typename Marker>
using FreeVarHolder = DisambiguationHolder<VET, Marker>;

namespace detail {
/**
 * @brief Get multiple free variables that can be combined in the same expression
 * 
 * @tparam MarkingType The unique type that is used to identify variables that are from the same expression
 * @param fnt the fixed number references to wrap in the variable expression
 */
template<typename MarkingType, FixedNumberType... FNT, std::size_t... indexes>
auto get_free_variables(FNT const & ... fnt, std::index_sequence<indexes...>) {
  return std::make_tuple(FreeVarHolder<Variable<FNT, indexes>, MarkingType>{fnt}...);
}
}

template<typename DisambiguationMarker = decltype([]{}), FixedNumberType... FNT>
auto FreeVariables(FNT const & ... fnt) {
  using indices = std::make_index_sequence<sizeof...(FNT)>;
  return detail::get_free_variables<DisambiguationMarker>(fnt..., indices{});
}

template<typename DisambiguationMarker = decltype([]{}), FixedNumberType FNT>
auto FreeVariable(FNT const & fnt) {
  return FreeVarHolder<Variable<FNT, 0>, DisambiguationMarker>{fnt};
}

namespace detail {
template <DisambiguationHType DT0, DisambiguationHType DT1>
static constexpr bool _are_holder_compatibles = std::is_same_v<typename DT0::disambiguation_t, typename DT1::disambiguation_t> || 
                                                std::is_same_v<typename DT0::disambiguation_t, detail::NoAmbiguities> || 
                                                std::is_same_v<typename DT1::disambiguation_t, detail::NoAmbiguities>;

template <DisambiguationHType DT0, DisambiguationHType DT1>
struct HolderUnionHelper {
  static_assert(_are_holder_compatibles<DT0, DT1>, "Holder union between incompatible types is impossible");
  using common_type = std::conditional_t<std::is_same_v<typename DT0::disambiguation_t, NoAmbiguities>, typename DT1::disambiguation_t, typename DT0::disambiguation_t>;
};                                              

template <DisambiguationHType DT0, DisambiguationHType DT1>
using _holder_union_t = typename HolderUnionHelper<DT0, DT1>::common_type;

template<UniqueExprType ET>
auto get_holder(ET const & expr) {
  return NonAmbiguousHolder<ET>{expr};
}

template<DisambiguationHType HT>
auto get_holder(HT const & holder) {
  return holder;
}
}

template<typename T>
concept ExprHoldType = requires (T val) { detail::get_holder(val); };

namespace detail {
template <template<ExpressionType, ExpressionType> class OperatorTType, ExprHoldType T1, ExprHoldType T2>
auto binary_op(T1 const & op1, T2 const & op2) {
  auto holder_1 = get_holder(op1);
  auto holder_2 = get_holder(op2);
  using h1_t = decltype(holder_1);
  using h2_t = decltype(holder_2);
  using expr_type = OperatorTType<typename h1_t::expression_t, typename h2_t::expression_t>;
  using ret_type = DisambiguationHolder<expr_type, detail::_holder_union_t<h1_t, h2_t>>; 
  return ret_type{expr_type{holder_1.expression, holder_2.expression}};
}

template <template<ExpressionType> class OperatorTType, ExprHoldType T1>
auto unary_op(T1 const & op) {
  auto holder = get_holder(op);
  using h_t = decltype(holder);
  using expr_type = OperatorTType<typename h_t::expression_t>;
  using return_t = DisambiguationHolder<expr_type, typename h_t::disambiguation_t>;
  return return_t{holder.expression};
}
}

} // namespace archgenlib

#endif // EXPRESSION_HPP
