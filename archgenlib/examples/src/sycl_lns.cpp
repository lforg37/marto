#include <cassert>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <ostream>
#include <string_view>
#include <type_traits>

#include "bitint_tools/type_helpers.hpp"
#include "fixedpoint/evaluator.hpp"
#include "hint.hpp"

#include "fixedpoint/constants.hpp"
#include "fixedpoint/expression.hpp"
#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"
#include "fixedpoint/literal.hpp"
#include "tools/static_math.hpp"

#include <iostream>

#ifndef GENERATING_HEADER
#include <sycl/sycl.hpp>
#endif

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

  storage_t exponent;
  int64_t value() const {
    return static_cast<uint64_t>(exponent.value());
  }
};

template <typename T> concept LNSType = T::is_lns_type;

template <LNSType OpT> constexpr OpT operator+(OpT const &op1, OpT const &op2) {
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
    auto expr_useful_bits =
        result_exp.template round_to<OpT::format.lsb_weight>();
    auto res = max_exp.modular_add(expr_useful_bits);
    return {.exponent = res};
  } else {
    return {.exponent = {0}};
  }
}

using lns_t = LNS<8, 8>;

int main() {
#ifdef GENERATING_HEADER
  lns_t a = {1 << 8};
  lns_t b = a + a;
  return 0;
#else
  constexpr std::size_t len = 3;
  sycl::buffer<lns_t, 1> a(len);
  sycl::queue Queue;

  {
    sycl::host_accessor a_a{a, sycl::write_only};
    a_a[0] = lns_t{1 << 8};
    a_a[1] = lns_t{1 << 8};
  }

  Queue.submit([&](sycl::handler &cgh) {
    sycl::accessor a_a{a, cgh, sycl::read_write};
    cgh.single_task<class FirstKernel>([=] {
      a_a[2] = a_a[0] + a_a[1];
    });
  });

  {
    sycl::host_accessor a_a{a, sycl::read_only};
    std::cout << a_a[2].value() << std::endl;
    assert(a_a[2].value() == 1 << 9);
  }
#endif
}
