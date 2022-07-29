#ifndef FIXEDPOINT_PLUGIN_GEN_EVALUATOR_HPP
#define FIXEDPOINT_PLUGIN_GEN_EVALUATOR_HPP

#ifndef ARCHGEN_USE_MLIR_PLUGIN
#error "Should be used with the MLIR plugin"
#endif

#include "operators.hpp"

#define ATTR_GEN_AS_MLIR __attribute__((annotate("archgen_mlir_emit_as_mlir")))

namespace archgenlib {

/// Never has a specialization header
constexpr bool has_specialization_header = false;

namespace detail {

struct ToBeFolded {};

template <typename T = ToBeFolded>
__attribute__((annotate("archgen_mlir_generic_op"))) T generic_op(...);

template <typename ET> struct evaluatorImpl {};

#define UNARY_OPERATOR(NAME, BACKEND_NAME, FRONTEND_NAME)                      \
  template <typename InnerET>                                                  \
  struct evaluatorImpl<UnaryOp<InnerET, OperationType<OperationKind::NAME>>> { \
    static ATTR_GEN_AS_MLIR ToBeFolded evaluate() {                            \
      return generic_op<ToBeFolded>(#BACKEND_NAME,                             \
                                    evaluatorImpl<InnerET>::evaluate());       \
    }                                                                          \
  };

#define BINARY_OPERATOR(NAME, BACKEND_NAME, FRONTEND_NAME)                     \
  template <typename LeftET, typename RightET>                                 \
  struct evaluatorImpl<                                                        \
      BinaryOp<LeftET, RightET, OperationType<OperationKind::NAME>>> {         \
    static ATTR_GEN_AS_MLIR ToBeFolded evaluate() {                            \
      return generic_op<ToBeFolded>(#BACKEND_NAME,                             \
                                    evaluatorImpl<LeftET>::evaluate(),         \
                                    evaluatorImpl<RightET>::evaluate());       \
    }                                                                          \
  };

#define NULLNARY_OPERATOR(NAME, BACKEND_NAME, FRONTEND_NAME)                   \
  template <> struct evaluatorImpl<NullaryOp<OperationType<NAME##Op>>> {       \
    static ATTR_GEN_AS_MLIR ToBeFolded evaluate() {                            \
      return generic_op<ToBeFolded>(#BACKEND_NAME);                            \
    }                                                                          \
  };


#include "operators.def"

template <typename NumTy, std::size_t ID>
struct evaluatorImpl<::archgenlib::Variable<NumTy, ID>> {
  static ATTR_GEN_AS_MLIR ToBeFolded evaluate() {
    return generic_op<ToBeFolded>("free_variable", generic_op<NumTy>("parameter", ID));
  }
};

template <typename FixedConstTy> struct evaluatorImpl<::archgenlib::Constant<FixedConstTy>> {
  static ATTR_GEN_AS_MLIR ToBeFolded evaluate() {
    return generic_op<ToBeFolded>("constant", ::archgenlib::FixedNumber<typename FixedConstTy::dimension_t>{FixedConstTy::value});
  }
};

template <typename T, typename ET, typename... Ts> ATTR_GEN_AS_MLIR T evaluateImpl(Ts... ts) {
  return detail::generic_op<T>("evaluate",
                               detail::evaluatorImpl<ET>::evaluate());
}

template <typename T, typename ET, typename... Ts> T evaluate(std::tuple<Ts...> var) {
  return std::apply([&](auto... val){
    return evaluateImpl<T, ET>(val...);
  }, var);
}

}
}

#undef ATTR_GEN_AS_MLIR
#endif
