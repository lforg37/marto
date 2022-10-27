#ifndef FIXEDPOINT_PLUGIN_GEN_EVALUATOR_HPP
#define FIXEDPOINT_PLUGIN_GEN_EVALUATOR_HPP

#ifndef ARCHGEN_USE_MLIR_PLUGIN
#error "Should be used with the MLIR plugin"
#endif

#include "operators.hpp"

#ifdef __SYCL_DEVICE_ONLY__
#define ARCHGEN_MLIR_EXTERNAL SYCL_EXTERNAL
#else
#define ARCHGEN_MLIR_EXTERNAL
#endif

#define ARCHGEN_MLIR_ATTR(STR) __attribute__((annotate("archgen_mlir_" #STR)))

namespace archgenlib {

/// Only one compilation always with implementation
constexpr bool has_implementation = true;

namespace detail {

struct ToBeFolded {};

template <typename T = ToBeFolded, typename... Ts>
ARCHGEN_MLIR_EXTERNAL
ARCHGEN_MLIR_ATTR(generic_op)
ARCHGEN_MLIR_ATTR(emit_as_mlir) T generic_op(Ts...);

template <typename ET> struct evaluatorImpl {};

#define UNARY_OPERATOR(NAME, BACKEND_NAME, FRONTEND_NAME)                      \
  template <typename InnerET>                                                  \
  struct evaluatorImpl<UnaryOp<InnerET, OperationType<OperationKind::NAME>>> { \
    static ARCHGEN_MLIR_ATTR(emit_as_mlir) ToBeFolded evaluate() {             \
      return generic_op<ToBeFolded>(#BACKEND_NAME,                             \
                                    evaluatorImpl<InnerET>::evaluate());       \
    }                                                                          \
  };

#define BINARY_OPERATOR(NAME, BACKEND_NAME, FRONTEND_NAME)                     \
  template <typename LeftET, typename RightET>                                 \
  struct evaluatorImpl<                                                        \
      BinaryOp<LeftET, RightET, OperationType<OperationKind::NAME>>> {         \
    static ARCHGEN_MLIR_ATTR(emit_as_mlir) ToBeFolded evaluate() {             \
      return generic_op<ToBeFolded>(#BACKEND_NAME,                             \
                                    evaluatorImpl<LeftET>::evaluate(),         \
                                    evaluatorImpl<RightET>::evaluate());       \
    }                                                                          \
  };

#define CONSTANT_OPERATOR(NAME, BACKEND_NAME, FRONTEND_NAME)                   \
  template <>                                                                  \
  struct evaluatorImpl<NullaryOp<OperationType<OperationKind::NAME>>> {        \
    static ARCHGEN_MLIR_ATTR(emit_as_mlir) ToBeFolded evaluate() {             \
      return generic_op<ToBeFolded>(#BACKEND_NAME);                            \
    }                                                                          \
  };

#include "operators.def"

template <typename NumTy, std::size_t ID>
struct evaluatorImpl<::archgenlib::Variable<NumTy, ID>> {
  static ARCHGEN_MLIR_ATTR(emit_as_mlir) ToBeFolded evaluate() {
    return generic_op<ToBeFolded>("variable",
                                  generic_op<NumTy>("parameter", ID));
  }
};

template <typename FixedConstTy>
struct evaluatorImpl<::archgenlib::Constant<FixedConstTy>> {
  static ARCHGEN_MLIR_ATTR(emit_as_mlir) ToBeFolded evaluate() {
    return generic_op<ToBeFolded>(
        "constant",
        ::archgenlib::FixedNumber<typename FixedConstTy::dimension_t>{
            FixedConstTy::value});
  }
};

/// everything here is passed by pointers such that the MLIR frontend doesn't
/// have to care about the LLVMIR frontend ABI decisions of how to pass
/// parameters and return types.
template <typename T, auto AM, typename ET, typename... Ts>
ARCHGEN_MLIR_ATTR(emit_as_mlir)
ARCHGEN_MLIR_ATTR(top_level) void evaluateImpl(T &res, Ts &...ts) {
  res = detail::generic_op<T>("evaluate",
                               detail::evaluatorImpl<ET>::evaluate(), AM);
}

template <typename T, auto AM, typename ET, typename... Ts>
T evaluate(std::tuple<Ts...> var) {
  T res;
  std::apply([&](auto... val) { return evaluateImpl<T, AM, ET>(res, val...); },
                    var);
  return res;
}

} // namespace detail
} // namespace archgenlib

#undef ARCHGEN_MLIR_ATTR
#endif
