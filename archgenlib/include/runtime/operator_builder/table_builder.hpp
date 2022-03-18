#ifndef RUNTIME_OPERATOR_BUILDER_TABLE_BUILDER_HPP
#define RUNTIME_OPERATOR_BUILDER_TABLE_BUILDER_HPP

#include <ostream>
#include <string_view>
#include <vector>

#include <gmpxx.h>

#include "runtime/expression_tree.hpp"

#include "operator.hpp"

namespace archgenlib {

namespace detail {
struct CPPTableInitialization : public detail::CPPExpr {
  std::vector<mpz_class> const values;
  FPDimRTRepr output_dim;
  virtual void to_stream(std::ostream &);
  CPPTableInitialization(std::vector<mpz_class> const &val,
                         FPDimRTRepr output_dim)
      : values{val}, output_dim{output_dim} {}
};
} // namespace detail

struct TableBuilder {
  FPDimRTRepr const &input_dim;
  FPDimRTRepr const &output_dim;
  std::vector<mpz_class> const &values;
  CPPOperator build_table();
  TableBuilder(FPDimRTRepr const &in_repr, FPDimRTRepr const &out_repr,
               std::vector<mpz_class> const &val)
      : input_dim{in_repr}, output_dim(out_repr), values{val} {}
};
} // namespace archgenlib

#endif