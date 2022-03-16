#ifndef RUNTIME_OPERATOR_BUILDER_TABLE_BUILDER_HPP
#define RUNTIME_OPERATOR_BUILDER_TABLE_BUILDER_HPP

#include <ostream>
#include <string_view>
#include <vector>

#include <gmpxx.h>

#include "runtime/expression_tree.hpp"

namespace archgenlib {

struct TableBuilder {
  std::vector<mpz_class> const values;
  FPDimRTRepr const input_dim;
  FPDimRTRepr const output_dim;

  inline void print_val(std::ostream &os, mpz_class const &val) const {
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
  }

  inline void stream_values(std::ostream &os) const {
    bool first = true;
    for (auto &val : values) {
      if (first) {
        first = false;
      } else {
        os << ",\n";
      }
      //os << "{ ";
      print_val(os, val);
      //os << " }";
    }
  }

  friend std::ostream &operator<<(std::ostream &os, TableBuilder const &tb) {
    os << "[] (auto e) {\n";
    os << "using ::hint::operator \"\"_ubi;\n"
       << "using ::hint::operator \"\"_sbi;\n"
       << "constexpr archgenlib::Table<" << tb.input_dim.width << ", "
       << tb.output_dim.toFPDimName() << "> values{";
    tb.stream_values(os);
    os << "};\n";
    os << "return values[e];\n";
    os << "}";
    return os;
  }
};
} // namespace archgenlib

#endif