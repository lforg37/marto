#ifndef RUNTIME_OPERATOR_BUILDER_AFFECTATION_HPP
#define RUNTIME_OPERATOR_BUILDER_AFFECTATION_HPP

#include <ostream>
#include <string_view>

namespace archgenlib {
template <typename SubExprT> struct Affectation {
  std::string_view const name;
  SubExprT const &sub_expr;
  auto get_name() const { return name; }
  friend std::ostream &operator<<(std::ostream &os, Affectation const &affect) {
    os << "auto " << affect.name << " = " << affect.sub_expr << ";\n";
    return os;
  }
  Affectation(std::string_view name, SubExprT const &sub_expr)
      : name{name}, sub_expr{sub_expr} {}
};
} // namespace archgenlib

#endif