#include <cassert>
#include <cstdint>
#include <ostream>
#include <string_view>
#include <type_traits>
#include <cmath>
#include <numbers>

#include "bitint_tools/type_helpers.hpp"
#include "hint.hpp"

#include "fixedpoint/expression.hpp"
#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"
#include "fixedpoint/constants.hpp"
#include "fixedpoint/literal.hpp"


#include <iostream>

template<typename FixedTy>
inline auto convert_to_double(FixedTy val) {
  auto val_d = static_cast<double>(val.value());
  val_d = std::ldexp(val_d, FixedTy::format_t::lsb_weight);
  return val_d;
}

template<int prec, typename FPTy, unsigned table_size>
struct wave_gen {
  using TableDim = archgenlib::FixedFormat<table_size - 1, 0, unsigned>;
  using TableFPTy = archgenlib::FixedNumber<TableDim>;
  using TableSizeCst = archgenlib::Constant<archgenlib::FixedConstant<
      archgenlib::FixedFormat<table_size, 0, unsigned>, 1 << table_size>>;
  FPTy val;
  FPTy step;
  wave_gen() : val{0}, step{0} {}
  wave_gen(FPTy phase, FPTy freq) : val{phase}, step{freq} {}
  auto lookup(TableFPTy v) {
    auto c =
        archgenlib::sin(0x2.p0_cst / TableSizeCst{} * archgenlib::pi *
                        archgenlib::Variable<TableFPTy>{v});
    return archgenlib::evaluate<prec>(c);
  }
  auto get_value(FPTy v) {
    /// if we dont need to interpolate simple lookup
    if constexpr (TableFPTy::format_t::lsb_weight <=
                  FPTy::format_t::lsb_weight) {
      return lookup(v);
    } else {
      TableFPTy high_bits = v.template extract<TableDim::msb_weight, 0>();
      auto low_bits = v.template extract<-1, FPTy::format_t::lsb_weight>();

      auto complement_low_bits = 0x1_fixed - low_bits;
      auto out1 = lookup(high_bits);
      auto out2 = lookup(high_bits.value() + 1);
      // std::cout << convert_to_double(out1) << " * "
      //           << convert_to_double(complement_low_bits) << " + "
      //           << convert_to_double(out2) << " * "
      //           << convert_to_double(low_bits) << " ";
      auto res = out1 * complement_low_bits + out2 * low_bits;
      return res;
    }
  }
  auto get_next() {
    auto res = get_value(val);
    val = {val.value() + step.value()};
    return res;
  }
};

#ifdef TARGET_VITIS
using fpdim_t = archgenlib::FixedFormat<9, -2, unsigned>;
using fpnum_t = archgenlib::FixedNumber<fpdim_t>;

using mul_t = archgenlib::FixedNumber<archgenlib::FixedFormat<-1, -8, unsigned>>;

template<int prec, typename FPTy, unsigned table_size, unsigned osc_log2 , unsigned max_freq>
struct additive_synt {
  static constexpr auto osc_count = 1 << osc_log2;
  std::array<wave_gen<prec, FPTy, table_size>, osc_count> oscs;
  additive_synt(int i) {
    for (int i = 0; i < oscs.size(); i++) {
      auto freq = (max_freq * (i + 1)) / osc_count;
      oscs[i] = wave_gen<prec, FPTy, table_size>(i * freq, freq);
    }
  }
  auto get_next(std::array<mul_t, 256> coef) {
    using rs_t = decltype(oscs[0].get_next());
    constexpr auto res_fmt = mul_t::format * rs_t::format;
    using acc_fmt = archgenlib::FixedFormat<res_fmt.msb_weight + osc_log2, res_fmt.lsb_weight, typename decltype(res_fmt)::sign_t>;
    using acc_t = archgenlib::FixedNumber<acc_fmt>;
    acc_t acc{0};
    for (int i = 0; i < oscs.size(); i++)
      acc += oscs[i].get_next() * coef[i];
    return acc;
  }
};

__VITIS_KERNEL auto test(int i, std::array<mul_t, 256> coef) {
  additive_synt<-5, archgenlib::FixedNumber<archgenlib::FixedFormat<9, -2, unsigned>>, 10, 8, 512> synt(i);
  return synt.get_next(coef);
}
#else

template<typename FixedTy>
bool compare(FixedTy val, double ref, int prec) {
  double to_double = convert_to_double(val);
  double interval = std::ldexp(double{1}, prec);
  bool is_in_interval = std::abs(to_double - ref) < interval;
  // std::cout << "compare: " << to_double << " ~= " << ref << " at " << interval << std::endl;
  if (!is_in_interval) {
    return false;
  }
  return true;
}

template<typename Format>
void check_freq(int offset, int f) {
  std::cout << "testing: " << __PRETTY_FUNCTION__ << " offset=" << offset << " freq=" << f << std::endl;

  archgenlib::FixedNumber<Format> freq = {f};
  wave_gen<-5, decltype(freq), 10> wave(offset, freq);
  int table_size = 1 << 10;
  double point = std::ldexp(double{1}, Format::lsb_weight) * offset;
  double step = std::ldexp(double{1}, Format::lsb_weight) * f;

  for (int i = 0; i < table_size; i++, point += step) {
    double ref = std::sin(2.0 * std::numbers::pi * point / (table_size));

    // std::cout << i << " ";
    if (!compare(wave.get_next(), ref, -5)) {
      assert(false);
    }
  }
}

int main() {
  if (!archgenlib::has_specialization_header)
    exit(0);
  check_freq<archgenlib::FixedFormat<9, 0, unsigned>>(0, 1);
  check_freq<archgenlib::FixedFormat<9, -2, unsigned>>(0, 1);
  check_freq<archgenlib::FixedFormat<9, 0, unsigned>>(0, 5);
  check_freq<archgenlib::FixedFormat<9, -2, unsigned>>(0, 5);
  check_freq<archgenlib::FixedFormat<9, 0, unsigned>>(45, 5);
  check_freq<archgenlib::FixedFormat<9, -2, unsigned>>(45, 5);
}
#endif
