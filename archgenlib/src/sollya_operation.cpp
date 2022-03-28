#include <cassert>

#include <gmpxx.h>
#include <mpfr.h>
#include <sollya.h>
#include <variant>

#include "runtime/expression_tree.hpp"
#include "fixedpoint/operations.hpp"
#include "runtime/sollya_operation.hpp"

namespace {
using namespace archgenlib;

sollya_obj_t build_sollya_node(ConstantLeafRTNode const & val);
sollya_obj_t build_sollya_node(VariableLeafRTNode const & val);
sollya_obj_t build_sollya_node(BinaryOpRTNode const & val);
sollya_obj_t build_sollya_node(UnaryOpRTNode const & val);
sollya_obj_t build_sollya_node(NullaryOpLeafRTNode const & val);

auto sollya_node_builder_visitor = []<typename T>(const T & val) {
  return build_sollya_node(val);
};

using sollya_binary_op_t = sollya_obj_t(sollya_obj_t, sollya_obj_t);
using sollya_unary_op_t = sollya_obj_t(sollya_obj_t);
using sollya_nullary_op_t = sollya_obj_t();

sollya_binary_op_t* get_bin_op(OperationKind op_kind) {
  assert(detail::isBinaryOpKind(op_kind));
  switch (op_kind) {
    case OperationKind::ADD:
    return sollya_lib_build_function_add;
    case OperationKind::SUB:
    return sollya_lib_build_function_sub;
    case OperationKind::MUL:
    return sollya_lib_build_function_mul;
    case OperationKind::DIV:
    return sollya_lib_build_function_div;
    case OperationKind::POW:
    return sollya_lib_build_function_pow;
    default:
    __builtin_unreachable();
  }
}

sollya_unary_op_t* get_unary_op(OperationKind op_kind) {
  assert(detail::isUnaryOpKind(op_kind));
  switch (op_kind) {
    case OperationKind::NEG:
      return sollya_lib_build_function_neg;
    case OperationKind::SIN:
      return sollya_lib_build_function_sin;
    case OperationKind::LOG2:
      return sollya_lib_build_function_log2;
    case OperationKind::LOG:
      return sollya_lib_build_function_log;
    case OperationKind::ABS:
      return sollya_lib_build_function_abs;
    default:
    __builtin_unreachable();
  }
}

sollya_nullary_op_t* get_nullary_op(OperationKind op_kind) {
  assert(detail::isNullaryOpKind(op_kind));
  switch (op_kind) {
    case OperationKind::PI:
      return sollya_lib_build_function_pi;
    default:
    __builtin_unreachable();
  }
}

sollya_obj_t
build_sollya_node(ConstantLeafRTNode const & val) {
  auto width = val->dim.width;
  mpz_class constant_repr{val->constant_representation};
  mpfr_t value;
  mpfr_init2(value, width);
  mpfr_set_z_2exp(value, constant_repr.get_mpz_t(), val->dim.lsb_weight,
                  MPFR_RNDN);
  auto ret = sollya_lib_constant(value);
  mpfr_clear(value);
  return ret;
}

sollya_obj_t build_sollya_node(VariableLeafRTNode const &) {
  return sollya_lib_build_function_free_variable();
}

sollya_obj_t build_sollya_node(BinaryOpRTNode const & binop) {
  auto left = std::visit(sollya_node_builder_visitor, binop->left);
  auto right = std::visit(sollya_node_builder_visitor, binop->right);
  auto sollya_func = get_bin_op(binop->operation_kind);
  return sollya_func(left, right);
}

sollya_obj_t build_sollya_node(UnaryOpRTNode const & op_node) {
  auto op = std::visit(sollya_node_builder_visitor, op_node->operand);
  auto sollya_func = get_unary_op(op_node->operation_kind);
  return sollya_func(op);
}

sollya_obj_t build_sollya_node(NullaryOpLeafRTNode const & op_node) {
  auto sollya_func = get_nullary_op(op_node->operation_kind);
  return sollya_func();
}

} // namespace
namespace archgenlib {
sollya_obj_t
sollyaFunctionFromNode(RTNode const &op_node) {
  return std::visit(sollya_node_builder_visitor, op_node);
}
} // namespace archgenlib
