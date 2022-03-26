#ifndef RUNTIME_OPERATOR_BUILDER_OPERATOR_HPP
#define RUNTIME_OPERATOR_BUILDER_OPERATOR_HPP

#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "runtime/expression_tree.hpp"

namespace archgenlib {

namespace detail {

struct CPPExpr;
using cpp_expr_ptr = std::shared_ptr<CPPExpr>;

class TypeSpecifier;

template<typename T>
concept ExprPtrT = std::is_same_v<std::remove_cvref_t<T>, cpp_expr_ptr>;

struct CPPExpr : std::enable_shared_from_this<CPPExpr> {
  virtual void to_stream(std::ostream &) = 0;
  friend std::ostream &operator<<(std::ostream &os, CPPExpr &expr) {
    expr.to_stream(os);
    return os;
  }
  cpp_expr_ptr assign_to(std::string_view name, TypeSpecifier ts);
  cpp_expr_ptr assign_to(std::string_view name);

  template<ExprPtrT... ArgTypes>
  cpp_expr_ptr call(std::string_view member_name, ArgTypes... args) {
    auto method = get(member_name);
    return (*method)(args...);
  }
  cpp_expr_ptr delegate_call_op(std::vector<cpp_expr_ptr> args);
  cpp_expr_ptr get(std::string_view prop_name);
  cpp_expr_ptr follow(std::string_view prop_name);
  cpp_expr_ptr operator[](cpp_expr_ptr sub_expr);
  template<ExprPtrT... ArgTypes>
  cpp_expr_ptr operator()(ArgTypes... args) {
    return delegate_call_op({args...});
  }
  cpp_expr_ptr to_return();
  cpp_expr_ptr use_to_build(std::string_view type_name);
};

struct TypeSpecifier {
  bool is_constexpr = false;
  bool is_static = false;
  std::optional<std::string> type_name;
  friend std::ostream& operator<<(std::ostream& os, TypeSpecifier const &);
};

struct Access : public CPPExpr {
  cpp_expr_ptr accessed;
  cpp_expr_ptr index;
  virtual void to_stream(std::ostream &os);
  Access(cpp_expr_ptr accessed, cpp_expr_ptr index)
      : accessed{accessed}, index{index} {}
};

struct Add : public CPPExpr {
  cpp_expr_ptr op0;
  cpp_expr_ptr op1;
  virtual void to_stream(std::ostream &os);
  Add(cpp_expr_ptr op0, cpp_expr_ptr op1)
      : op0{op0}, op1{op1} {}
};

struct Assignment : public CPPExpr {
  bool is_declaration;
  std::string name;
  cpp_expr_ptr expr;
  TypeSpecifier type_specifier;
  virtual void to_stream(std::ostream &os);
  Assignment(std::string_view name, cpp_expr_ptr expr, TypeSpecifier ts = {})
      : name{name}, expr{expr}, is_declaration{true}, type_specifier{ts} {}
};

struct Call : public CPPExpr {
  cpp_expr_ptr called;
  std::vector<cpp_expr_ptr> args;
  virtual void to_stream(std::ostream &os);
  Call(cpp_expr_ptr called, std::vector<cpp_expr_ptr> args= {})
      : called{called}, args{args} {}

};

struct Construction : public CPPExpr {
  std::string type_desc;
  cpp_expr_ptr expr;
  virtual void to_stream(std::ostream &os);
  Construction(std::string_view type_name, cpp_expr_ptr expr)
      : type_desc{type_name}, expr{expr} {}
};

struct InputDesc : public CPPExpr {
  bool is_declaration;
  std::string name;
  std::optional<FixedFormatRTRepr> dimension;
  virtual void to_stream(std::ostream &os);
  InputDesc(std::string_view name):name{name}, is_declaration(true){}
  InputDesc(std::string_view name, FixedFormatRTRepr dim):name{name},dimension{dim},is_declaration{true}{}
};

struct MemberCall : public CPPExpr {
  cpp_expr_ptr expression;
  bool should_dereference;
  std::string property_name;
  virtual void to_stream(std::ostream &);
  MemberCall(cpp_expr_ptr expr, bool call_is_indirect, std::string_view pname)
      : expression{expr}, should_dereference{call_is_indirect}, property_name{
                                                                    pname} {}
};

struct NamedExpression : public CPPExpr {
  std::string name;
  virtual void to_stream(std::ostream& os);
  NamedExpression(std::string_view name):name{name}{}
};

struct Return : public CPPExpr {
  cpp_expr_ptr expression;
  virtual void to_stream(std::ostream &);
  Return(cpp_expr_ptr expr) : expression{expr} {}
};

} // namespace detail

struct CPPOperator : public detail::CPPExpr {
  // TODO find something cleaner to ensure only InputDesc are here
  std::vector<detail::cpp_expr_ptr> inputs;

  std::optional<FixedFormatRTRepr> output_dim;
  std::vector<detail::cpp_expr_ptr> instructions;
  virtual void to_stream(std::ostream &);
};

detail::cpp_expr_ptr operator+(detail::cpp_expr_ptr op0, detail::cpp_expr_ptr op1); 

} // namespace archgenlib
#endif
