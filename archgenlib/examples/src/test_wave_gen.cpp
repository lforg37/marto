#include <cstdint>
#include <string_view>
#include <type_traits>
#include <numbers>
#include <bit>

#include "bitint_tools/type_helpers.hpp"
#include "hint.hpp"

#include "fixedpoint/expression.hpp"
#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"
#include "fixedpoint/constants.hpp"
#include "fixedpoint/literal.hpp"

template <bool fail, auto...> struct print_assert {
  static_assert(fail, "print");
};

template <auto... vals> struct print : print_assert<false, vals...> {};

using archgenlib::bitweight_t;  

constexpr bitweight_t outprec = -10;

using fpdim_t = archgenlib::FixedFormat<-1, -14, unsigned>;
using fpnum_t = archgenlib::FixedNumber<fpdim_t>;

using storage_t =
    hint::detail::bitint_base_t<fpdim_t::is_signed, fpdim_t::width>;

#ifndef __VITIS_KERNEL
#warning "unsupported compiler"
#define __VITIS_KERNEL
#endif

template<int prec, typename FPTy, unsigned table_size>
struct wave_gen {
  using TableFPTy = archgenlib::FixedNumber<
      archgenlib::FixedFormat<-1, -static_cast<int>(table_size), unsigned>>;
  FPTy val;
  FPTy step;
  wave_gen(FPTy f) : val{0}, step{0} {
    set_freq(f);
  }
  void set_freq(FPTy freq) {
    step = freq;
  }
  auto lookup(TableFPTy v) {
    auto c = archgenlib::sin(0x2.p0_cst * archgenlib::pi *
                             archgenlib::Variable<TableFPTy>{v});
    return archgenlib::evaluate<prec>(c);
  }
  auto get_value(FPTy v) {
    /// if we dont need to interpolate simple lookup
    if constexpr (TableFPTy::format_t::lsb_weight <=
                  FPTy::format_t::lsb_weight) {
      return lookup(v);
    } else {
      auto bits = v.as_hint();
      using bit_wrapper = decltype(bits);
      TableFPTy high_bits{
          bits.template slice<bit_wrapper::width - 1,
                              bit_wrapper::width - table_size>()};
      archgenlib::FixedNumber<archgenlib::FixedFormat<
          -1, -((int)bit_wrapper::width - (int)table_size), unsigned>>
          low_bits{bits.template slice<bit_wrapper::width - table_size - 1, 0>()};

      auto complement_low_bits = 0x1.p0_fixed - low_bits;
      auto out1 = lookup(high_bits);
      auto out2 = lookup(high_bits.value() + 1);
      return out1 * low_bits + out2 * complement_low_bits;
    }
  }
  auto get_next() {
    auto res = get_value(val);
    val = {val.value() + step.value()};
    return res;
  }
};

__VITIS_KERNEL auto test(int i) {
  wave_gen<-5, fpnum_t, 10> wave(1 << 5);
  return wave.get_next();
}
