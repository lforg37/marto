#include <iostream>
#include <memory>
#include <ostream>
#include <sstream>
#include <string_view>

#include "fixedpoint/fixedpoint.hpp"
#include "fixedpoint/operators/table.hpp"
#include "runtime/expression_tree.hpp"
#include "runtime/operator_builder/multipartite_operator.hpp"
#include "runtime/operator_builder/operator.hpp"
#include "runtime/operator_builder/table_builder.hpp"
#include "runtime/operator_builder/value_getter.hpp"
#include "runtime/output_formatter.hpp"

namespace archgenlib {
namespace detail {
cpp_expr_ptr CPPExpr::assign_to(std::string_view name) {
  return cpp_expr_ptr{new Assignment(name, shared_from_this())};
}

cpp_expr_ptr CPPExpr::assign_to(std::string_view name, TypeSpecifier ts) {
  return cpp_expr_ptr{new Assignment(name, shared_from_this(), ts)};
}

cpp_expr_ptr CPPExpr::get(std::string_view prop_name) {
  return cpp_expr_ptr{new MemberCall(shared_from_this(), false, prop_name)};
}

cpp_expr_ptr CPPExpr::follow(std::string_view prop) {
  return cpp_expr_ptr{new MemberCall(shared_from_this(), true, prop)};
}

cpp_expr_ptr CPPExpr::operator[](cpp_expr_ptr sub_expr) {
  return cpp_expr_ptr{new Access{shared_from_this(), sub_expr}};
}

cpp_expr_ptr CPPExpr::delegate_call_op(std::vector<cpp_expr_ptr> sub_expr) {
  return cpp_expr_ptr{new Call{shared_from_this(), sub_expr}};
}

cpp_expr_ptr CPPExpr::to_return() {
  return cpp_expr_ptr{new Return{shared_from_this()}};
}

cpp_expr_ptr CPPExpr::use_to_build(std::string_view type_desc) {
  return cpp_expr_ptr{new Construction(type_desc, shared_from_this())};
}

std::ostream &operator<<(std::ostream &os, TypeSpecifier const &ts) {
  if (ts.is_constexpr) {
    os << "constexpr ";
  }
  if (ts.type_name.has_value()) {
    os << *ts.type_name;
  } else {
    os << "auto";
  }
  return os;
}

void Access::to_stream(std::ostream &os) {
  os << *accessed << "[" << *index << "]";
}

void Add::to_stream(std::ostream &os) { os << *op0 << "+" << *op1; }

void Assignment::to_stream(std::ostream &os) {
  if (is_declaration) {

    os << type_specifier << " " << name << " = " << *expr << ";\n";
    is_declaration = false;
  } else {
    os << name;
  }
}

void Call::to_stream(std::ostream &os) {
  os << *called << "(";
  bool is_first = true;
  for (auto &arg : args) {
    if (is_first) {
      is_first = false;
    } else {
      os << ", ";
    }
    os << *arg;
  }
  os << ")";
}

void Construction::to_stream(std::ostream &os) {
  os << type_desc << "{" << *expr << "}";
}

void CPPTableInitialization::to_stream(std::ostream &os) {
  os << "{\n";
  bool first = true;

  auto stream_val = [&os, this](mpz_class const &val) {
    mpz_class to_encode;
    mpz_class mask{1};
    mask <<= output_dim.width;
    if (val < 0) {
      assert(output_dim.is_signed);
      to_encode = mask + val;
    } else {
      to_encode = val;
    }
    os << "0b";
    for (mask >>= 1; mask > 0; mask >>= 1) {
      if ((to_encode & mask) > 0) {
        os << "1";
      } else {
        os << "0";
      }
    }
    if (output_dim.is_signed) {
      os << "_sbi";
    } else {
      os << "_ubi";
    }
  };

  for (auto &val : values) {
    if (first) {
      first = false;
    } else {
      os << ",\n";
    }
    stream_val(val);
  }
  os << "}";
}

void InputDesc::to_stream(std::ostream &os) {
  if (is_declaration) {
    if (dimension.has_value()) {
      os << dimension->toFPNumName() << " ";
    } else {
      os << "auto ";
    }
    is_declaration = false;
  }
  os << name;
}

void MemberCall::to_stream(std::ostream &os) {
  os << *expression;
  if (should_dereference) {
    os << "->";
  } else {
    os << ".";
  }
  os << property_name;
}

void NamedExpression::to_stream(std::ostream &os) {
  os << name;
}

void Return::to_stream(std::ostream &os) {
  os << "return " << *expression << ";\n";
}

} // namespace detail

void CPPOperator::to_stream(std::ostream &os) {
  os << "[]";
  if (!inputs.empty()) {
    os << "(";
    bool first = true;
    for (auto &input_ptr : inputs) {
      if (!first)
        os << ", ";
      input_ptr->to_stream(os);
    }
    os << ") ";
  }
  if (output_dim.has_value()) {
    os << "->" << output_dim->toFixedFormatName() << " ";
  }
  os << "{\n";
  for (auto &instr_ptr : instructions) {
    instr_ptr->to_stream(os);
  }
  os << "}";
}

detail::cpp_expr_ptr get_path(ExpressionRTRepr::path_t const &path,
                              detail::cpp_expr_ptr main_expr) {
  for (auto iter = path.rbegin(); iter != path.rend(); iter++) {
    std::string_view direction;
    switch (*iter) {
    case ExpressionRTRepr::DirectionInTree::LEFT:
      direction = "left";
      break;
    case ExpressionRTRepr::DirectionInTree::RIGHT:
      direction = "right";
      break;
    case ExpressionRTRepr::DirectionInTree::DIRECT:
      direction = "op";
      break;
    }
    main_expr = main_expr->get(direction);
  }
  main_expr = main_expr->get("value");
  return main_expr;
}

CPPOperator TableBuilder::build_table() {
  CPPOperator ret{};
  std::stringstream s{};

  detail::cpp_expr_ptr input{new detail::InputDesc{"point", input_dim}};
  ret.inputs.emplace_back(input);
  s << "archgenlib::Table<" << input_dim.width << ", "
    << output_dim.toFixedFormatName() << ">";
  detail::cpp_expr_ptr tab_init_list{
      new detail::CPPTableInitialization(values, output_dim)};
  auto tab_val = tab_init_list->use_to_build(s.str())->assign_to(
      "values", detail::TypeSpecifier{.is_constexpr = true});
  ret.instructions.emplace_back(tab_val);
  auto tab_key = input->get("value()")->assign_to("to_num_key");
  ret.instructions.emplace_back(tab_key);
  auto res = (*tab_val)[tab_key]->to_return();
  ret.instructions.emplace_back(res);
  return ret;
}

CPPOperator MultipartiteOperator::get_operator() {
  assert(mpf.best_config.has_value());
  assert(mpf.check_best_config(-38));
  CPPOperator ret{};
  detail::cpp_expr_ptr input{
      new detail::InputDesc{"point", mpf.function.input_format}};
  ret.inputs.emplace_back(input);

  auto slicer = [](detail::cpp_expr_ptr val, unsigned int MSBRel,
                   unsigned int LSBRel) {
    std::stringstream ss{};
    ss << "template slice<" << MSBRel << ", " << LSBRel << ">";
    return val->call(ss.str());
  };

  FixedFormatRTRepr output_dim{mpf.function.msb_output, mpf.lsb_out, mpf.function.signed_output};

  // TIV
  auto &config = mpf.best_config.value();
  FixedFormatRTRepr table_output_dim{mpf.function.msb_output, mpf.lsb_out - 2, // Two guard bits
                               mpf.function.signed_output};
  FixedFormatRTRepr tiv_input_dim{static_cast<bitweight_t>(config.alpha) - 1, 0,
                            false};
  auto tiv_builder = TableBuilder{tiv_input_dim, table_output_dim,
                                  config.initial_values_table};
  detail::cpp_expr_ptr tiv_op{new CPPOperator{tiv_builder.build_table()}};
  auto tiv_func = tiv_op->assign_to("initial_value");
  ret.instructions.emplace_back(tiv_func);
  auto in_as_hint = input->call("as_hint")->assign_to("hint_input");
  ret.instructions.emplace_back(in_as_hint);
  auto alpha_slice = slicer(in_as_hint, mpf.function.input_format.width - 1,
                            mpf.function.input_format.width - config.alpha);
  auto tiv_val =
      (*tiv_func)(alpha_slice->call("unravel"))->assign_to("tiv_val");
  ret.instructions.emplace_back(tiv_val);

  // TO - TODO : replace when going to multipartite;
  auto to_key_high_bits = slicer(
      in_as_hint, mpf.function.input_format.width - 1,
      mpf.function.input_format.width - config.offset_configs[0].gamma);
  auto to_key_low_bits =
      slicer(in_as_hint, config.offset_configs[0].beta - 1, 0);
  auto to_key = to_key_high_bits->call("concatenate", to_key_low_bits)
                    ->call("unravel")
                    ->assign_to("to_key");
  ret.instructions.emplace_back(to_key);
  auto &to_config = config.offset_configs[0];
  FixedFormatRTRepr to_input_dim{
      static_cast<bitweight_t>(to_config.gamma + to_config.beta) - 1, 0, false};
  auto to_builder =
      TableBuilder{to_input_dim, table_output_dim, config.offset_tables[0]};
  detail::cpp_expr_ptr to_op{new CPPOperator{to_builder.build_table()}};
  auto to_func = to_op->assign_to("offsets");
  ret.instructions.emplace_back(to_func);

  auto to_val = (*to_func)(to_key)->assign_to("to_value");
  ret.instructions.emplace_back(to_val);

  auto res_val_scaled = to_val->call("as_hint")->call("as_unsigned")->call(
      "modularAdd",
      tiv_val->call("as_hint")->call(
          "modularAdd",
          detail::cpp_expr_ptr{new detail::NamedExpression{"{1}"}}));
  auto sliced_res = slicer(res_val_scaled, table_output_dim.width - 1, 2); // Guard bits
  auto retval = sliced_res->call("unravel")->use_to_build(output_dim.toFPNumName())->to_return();
  ret.instructions.emplace_back(retval);
  return ret;
}

detail::cpp_expr_ptr operator+(detail::cpp_expr_ptr op0,
                               detail::cpp_expr_ptr op1) {
  return detail::cpp_expr_ptr{new detail::Add{op0, op1}};
}

}; // namespace archgenlib
