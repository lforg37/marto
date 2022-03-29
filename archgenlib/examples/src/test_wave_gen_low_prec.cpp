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

template <typename FixedTy> inline auto convert_to_double(FixedTy val) {
  auto val_d = static_cast<double>(val.value());
  val_d = std::ldexp(val_d, FixedTy::format_t::lsb_weight);
  return val_d;
}

template <int prec, typename FPTy, unsigned table_size, auto f> struct wave_gen {
  using TableDim = archgenlib::FixedFormat<table_size - 1, 0, unsigned>;
  using TableFPTy = archgenlib::FixedNumber<TableDim>;
  using TableSizeCst = archgenlib::Constant<archgenlib::FixedConstant<
      archgenlib::FixedFormat<table_size-1, table_size-1, unsigned>, 1>>;
  FPTy phase;
  static constexpr auto freq = FPTy{f};
  wave_gen(FPTy phase) : phase{phase} {
  }
  auto lookup(TableFPTy v) {
    auto table_size_cst = TableSizeCst{};
    auto var = archgenlib::Variable{v};
    auto expr = archgenlib::sin( archgenlib::pi  * archgenlib::Variable<TableFPTy>{v} / 
                                  table_size_cst);
    return archgenlib::evaluate<prec>(expr);
  }

  auto get_value(FPTy v) {
    /// if we dont need to interpolate simple lookup
    // static_assert(TableFPTy::format_t::lsb_weight <=
    //               FPTy::format_t::lsb_weight);
    if constexpr (TableFPTy::format_t::lsb_weight <=
                  FPTy::format_t::lsb_weight) {
      return lookup(v);
    } else {
      TableFPTy high_bits = v.template extract<TableDim::msb_weight, 0>();
      auto low_bits = v.template extract<-1, FPTy::format_t::lsb_weight>();
      auto complement_low_bits = (0x1_fixed - low_bits).template extract<0, FPTy::format_t::lsb_weight>();
      auto out1 = lookup(high_bits).modular_mult(complement_low_bits);
      auto out2 = lookup(high_bits.value() + 1).modular_mult(low_bits);
      // std::cout << convert_to_double(out1) << " * "
      //           << convert_to_double(complement_low_bits) << " + "
      //           << convert_to_double(out2) << " * "
      //           << convert_to_double(low_bits) << " " << std::endl;
     
      return out1.modular_add(out2);
      // auto res = out1 * complement_low_bits + out2 * low_bits;
      // std::cout << convert_to_double(res) << std::endl;
      // return res;
    }
  }
  auto get_next() {
    auto res = get_value(phase);
    phase += freq;
    return res;
  }
};

namespace detail {
template<unsigned start_idx, unsigned last_idx, unsigned int NBOsc, int prec, typename FPTy, unsigned table_size, auto max_frequency, typename Enable = void>
struct OscillatorBench;

template<unsigned start_idx, unsigned last_idx, unsigned int NBOsc, int prec, typename FPTy, unsigned table_size, auto max_frequency>
struct OscillatorBench<start_idx, last_idx, NBOsc, prec, FPTy, table_size, max_frequency, std::enable_if_t<start_idx == last_idx>>
{
  static constexpr auto freq = max_frequency * start_idx / NBOsc;
  using osc_t = wave_gen<prec, FPTy, table_size, freq>;
  osc_t oscillator{0};
  template<archgenlib::FixedNumberType T>
  auto result(const std::array<T, NBOsc>& coef) {
    constexpr auto guard_bits = std::min(hint::Static_Val<NBOsc>::_log2, T::width);
    auto sinval = oscillator.get_next();
    // std::cout << convert_to_double(sinval) << std::endl;
    auto prod = sinval.modular_mult(coef[start_idx - 1]);
    // using prod_t = decltype(prod);
    // constexpr auto low = std::max<archgenlib::bitweight_t>(prec-guard_bits, prod_t::format.lsb_weight); 
    // auto res = prod.template extract<prod_t::format.msb_weight, low>();
    return prod;
  }
};

template<unsigned start_idx, unsigned last_idx, unsigned int NBOsc, int prec, typename FPTy, unsigned table_size, auto max_frequency>
struct OscillatorBench<start_idx, last_idx, NBOsc, prec, FPTy, table_size, max_frequency, std::enable_if_t<start_idx < last_idx>>
{
  static constexpr unsigned mid_point = (last_idx + start_idx) / 2;
  using low_type = OscillatorBench<start_idx, mid_point, NBOsc, prec, FPTy, table_size, max_frequency>;
  using high_type = OscillatorBench<mid_point + 1, last_idx, NBOsc, prec, FPTy, table_size, max_frequency>;
  low_type low{};
  high_type high{};
  auto result(const auto& coef) {
    auto res_low = low.result(coef);
    auto res = res_low.modular_add(high.result(coef));\
    // std::cout << convert_to_double(res) << std::endl;
    return res;
  }
};

}


template<unsigned int NBOsc, int prec, typename FPTy, unsigned table_size, auto max_frequency>
struct AdditiveSynthesizer {
  using osc_bank_t = detail::OscillatorBench<1, NBOsc, NBOsc, prec, FPTy, table_size, max_frequency>;
  osc_bank_t osc_bank{};
  auto get_value(const auto& coef) {
    auto unrounded_res = osc_bank.result(coef);
    // std::cout << convert_to_double(unrounded_res) << std::endl;
    return unrounded_res;
    // return unrounded_res.template round_to<prec>();
  }
};

template<unsigned int NBOsc, int prec, typename FPTy, unsigned table_size, auto max_frequency>
AdditiveSynthesizer<NBOsc, prec, FPTy, table_size, max_frequency> additive_synth{};

// using fpdim_t = archgenlib::FixedFormat<9, -2, unsigned>;
// using fpnum_t = archgenlib::FixedNumber<fpdim_t>;

using mul_t =
    archgenlib::FixedNumber<archgenlib::FixedFormat<-1, -8, unsigned>>;

template<typename T, auto size>
__attribute((always_inline)) auto test2(std::array<T, size>& coef) {
  using fixe_t = archgenlib::FixedNumber<archgenlib::FixedFormat<7, 0, unsigned>>;
  auto& add_synth = additive_synth<size, -5, fixe_t, 8, 255>;
  auto res =  add_synth.get_value(coef);
  // std::cout << "test2:" << convert_to_double(res) << std::endl;
  return res;
}

#ifdef TARGET_VITIS
__VITIS_KERNEL auto test(std::array<mul_t, 256> coef) {
  static_assert(archgenlib::has_specialization_header);
  return test2(coef);
}
#else

template <typename FixedTy> bool compare(FixedTy val, double ref, int prec) {
  double to_double = convert_to_double(val);
  double interval = std::ldexp(double{1}, prec);
  bool is_in_interval = std::abs(to_double - ref) < interval;
  // std::cout << "compare: " << to_double << " ~= " << ref << " at " <<
  // interval << std::endl;
  if (!is_in_interval) {
    return false;
  }
  return true;
}

template <typename Format, int f> void check_freq(int offset) {
  std::cout << "testing: " << __PRETTY_FUNCTION__ << " offset=" << offset
            << " freq=" << f << std::endl;

  using fp_num = archgenlib::FixedNumber<Format>;
  wave_gen<-5, fp_num , 10, f> wave(offset);
  int table_size = 1 << 10;
  double point = std::ldexp(double{1}, Format::lsb_weight) * offset;
  double step = std::ldexp(double{1}, Format::lsb_weight) * f;

  for (int i = 0; i < table_size; i++, point += step) {
    double ref = std::sin(2.0 * std::numbers::pi * point / (table_size));

    // std::cout << i << " ";
    if (!compare(wave.get_from_time({i}), ref, -5)) {
      assert(false);
    }
  }
}

void plot() {
  std::array<mul_t, 2> coef = {};
  coef[0] = mul_t::get_from_value(0.5);
  coef[1] = mul_t::get_from_value(.5);
  for (int i = 0; i < (1 << 8); i += 1) {
    auto res = convert_to_double(test2(coef));
    std::cout << i << "," << res << std::endl;
  }
}

int main() {
  if (!archgenlib::has_specialization_header)
    exit(0);
  plot();
  // check_freq<archgenlib::FixedFormat<9, 0, unsigned>, 1>(0);
  // check_freq<archgenlib::FixedFormat<9, -2, unsigned>, 1>(0);
  // check_freq<archgenlib::FixedFormat<9, 0, unsigned>, 5>(0);
  // check_freq<archgenlib::FixedFormat<9, -2, unsigned>, 5>(0);
  // check_freq<archgenlib::FixedFormat<9, 0, unsigned>, 5>(45);
  // check_freq<archgenlib::FixedFormat<9, -2, unsigned>, 5>(45);
}
#endif
