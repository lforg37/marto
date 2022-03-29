#include <cassert>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <ostream>
#include <string_view>
#include <type_traits>

#include "bitint_tools/type_helpers.hpp"
#include "hint.hpp"

#include "fixedpoint/constants.hpp"
#include "fixedpoint/expression.hpp"
#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"
#include "fixedpoint/literal.hpp"
#include "tools/static_math.hpp"

#include <iostream>
using namespace archgenlib;

/**
 * @brief Represent a LNS format
 *
 * @tparam NIntBits
 * @tparam NFracBits
 */
template <unsigned int NIntBits, unsigned int NFracBits> struct LNS {
  static constexpr bool is_lns_type = true;
  static_assert(NIntBits != 0 || NFracBits != 0);
  using exp_storage_fmt =
      archgenlib::FixedFormat<NIntBits - 1, -static_cast<int>(NFracBits),
                              unsigned>;
  using storage_t = FixedNumber<exp_storage_fmt>;
  static constexpr exp_storage_fmt format{};

  bool is_zero;
  bool is_neg;
  storage_t exponent;

  static constexpr LNS zero() { return {.is_zero = true}; }
};

template <typename T> concept LNSType = T::is_lns_type;

template <LNSType OpT> constexpr OpT operator+(OpT const &op1, OpT const &op2) {
  if (op1.is_zero)
    return op2;
  if (op2.is_zero)
    return op1;
  assert(op1.is_neg == op2.is_neg);
  typename OpT::storage_t max_exp, min_exp;
  bool op1_is_max;
  if (op1.exponent > op2.exponent) {
    op1_is_max = true;
    max_exp = op1.exponent;
    min_exp = op2.exponent;
  } else {
    op1_is_max = false;
    max_exp = op2.exponent;
    min_exp = op1.exponent;
  }
  auto exp_diff = min_exp.modular_sub(max_exp);
  auto diff_var = Variable{exp_diff};
  auto expr = archgenlib::log2(0x1.p0_cst + pow(0x2.p0_cst, diff_var));

  constexpr auto rounding_bit_weight = OpT::format.lsb_weight - 1;
  using round_bit_format =
      FixedFormat<rounding_bit_weight, rounding_bit_weight, unsigned>;

  // Rounding in log space
  auto expr_rounded = expr + Constant<FixedConstant<round_bit_format, 1>>{};

  auto result_exp = evaluate<rounding_bit_weight - 1>(expr_rounded);
  // auto extract_res = result_exp.template round_to<OpT::format.lsb_weight>();
  if constexpr (has_specialization_header) {
    using expr_formt = typename decltype(result_exp)::format_t;
    auto expr_useful_bits = result_exp.template round_to<OpT::format.lsb_weight>();
    auto res = max_exp.modular_add(expr_useful_bits);
    return {.is_zero = false, .is_neg = op1.is_neg, .exponent = res};
  } else {
    return {.exponent = {0}};
  }
}

using lns_t = LNS<6, 6>;

#ifdef TARGET_VITIS
__VITIS_KERNEL auto test(lns_t a, lns_t b) {
  return a + b;
}

#else

int main() {
  lns_t a{.is_zero = false, .is_neg = false, .exponent{0b10000}};
  auto res = a + a;
  std::cout << (int)res.exponent.value() << std::endl;
}
#endif
