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
  FixedFormatRTRepr output_dim;
  virtual void to_stream(std::ostream &);
  CPPTableInitialization(std::vector<mpz_class> const &val,
                         FixedFormatRTRepr output_dim)
      : values{val}, output_dim{output_dim} {}
};
} // namespace detail

struct TableBuilder {
  FixedFormatRTRepr const &input_dim;
  FixedFormatRTRepr const &output_dim;
  std::vector<mpz_class> const &values;
  CPPOperator build_table();
  TableBuilder(FixedFormatRTRepr const &in_repr, FixedFormatRTRepr const &out_repr,
               std::vector<mpz_class> const &val)
      : input_dim{in_repr}, output_dim(out_repr), values{val} {}
};
} // namespace archgenlib

#endif