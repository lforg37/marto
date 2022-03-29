#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <ostream>
#include <string_view>
#include <type_traits>


#include <iostream>

template <typename FixedTy> inline auto convert_to_double(FixedTy val) {
  auto val_d = static_cast<double>(val.value());
  val_d = std::ldexp(val_d, FixedTy::format_t::lsb_weight);
  return val_d;
}

template <auto fNum, auto fDenom> struct wave_gen {
  float phase;
  static constexpr auto freq = static_cast<float>(fNum)/fDenom;
  wave_gen(float phase) : phase{phase} {}

  auto get_from_time(float time) {
    time *= 2*std::numbers::pi_v<float> * freq;
    time += phase;
    return std::sin(time);
  }
};

namespace detail {
template<unsigned start_idx, unsigned last_idx, unsigned int NBOsc, auto max_frequency, typename Enable = void>
struct OscillatorBench;

template<unsigned start_idx, unsigned last_idx, unsigned int NBOsc, auto max_frequency>
struct OscillatorBench<start_idx, last_idx, NBOsc, max_frequency, std::enable_if_t<start_idx == last_idx>>
{
  static constexpr auto freqNum = max_frequency * start_idx / NBOsc;
  using osc = wave_gen<freqNum, max_frequency>;

  auto result(float inval, const std::array<float, NBOsc>& coef) {
    return coef[start_idx - 1] * osc{0.}.get_from_time(inval);
  }
};

template<unsigned start_idx, unsigned last_idx, unsigned int NBOsc, auto max_frequency>
struct OscillatorBench<start_idx, last_idx, NBOsc, max_frequency, std::enable_if_t<start_idx < last_idx>>
{
  static constexpr unsigned mid_point = (last_idx + start_idx) / 2;
  using low_type = OscillatorBench<start_idx, mid_point, NBOsc, max_frequency>;
  using high_type = OscillatorBench<mid_point + 1, last_idx, NBOsc, max_frequency>;
  auto result(float inval, const auto& coef) {
    auto res_low = low_type{}.result(inval, coef);
    return res_low + high_type{}.result(inval, coef);
  }
};

}

template<unsigned int NBOsc, auto max_frequency>
struct AdditiveSynthesizer {
  using ul_type = detail::OscillatorBench<1, NBOsc, NBOsc, max_frequency>;
  auto get_value(float time, const auto& coef) {
    return ul_type{}.result(time);
  }
};

#ifdef TARGET_VITIS
using fpdim_t = archgenlib::FixedFormat<9, -2, unsigned>;
using fpnum_t = archgenlib::FixedNumber<fpdim_t>;

using mul_t =
    archgenlib::FixedNumber<archgenlib::FixedFormat<-1, -8, unsigned>>;

// template <int prec, typename FPTy, unsigned table_size, unsigned osc_log2,
//           unsigned max_freq>
// struct additive_synt {
//   static constexpr auto osc_count = 1 << osc_log2;
//   std::array<wave_gen<prec, FPTy, table_size>, osc_count> oscs;
//   additive_synt(int i) {
//     for (int i = 0; i < oscs.size(); i++) {
//       auto freq = (max_freq * (i + 1)) / osc_count;
//       oscs[i] = wave_gen<prec, FPTy, table_size>(i * freq, freq);
//     }
//   }
//   auto get_next(std::array<mul_t, 256> coef) {
//     using rs_t = decltype(oscs[0].get_next());
//     constexpr auto res_fmt = mul_t::format * rs_t::format;
//     using acc_fmt = archgenlib::FixedFormat<res_fmt.msb_weight + osc_log2,
//                                             res_fmt.lsb_weight,
//                                             typename decltype(res_fmt)::sign_t>;
//     using acc_t = archgenlib::FixedNumber<acc_fmt>;
//     acc_t acc{0};
//     for (int i = 0; i < oscs.size(); i++)
//       acc += oscs[i].get_next() * coef[i];
//     return acc;
//   }
// };

// unsigned int NBOsc, int prec, typename FPTy, unsigned table_size, auto max_frequency

__VITIS_KERNEL auto test(int i, std::array<mul_t, 256> coef) {
  using fixe_t = archgenlib::FixedNumber<archgenlib::FixedFormat<9, -2, unsigned>>;
  AdditiveSynthesizer<256, -5, fixe_t, 10, 4 << 2> synt;
  static_assert(archgenlib::has_specialization_header);
  // additive_synt<
  //     -5, archgenlib::FixedNumber<archgenlib::FixedFormat<9, -2, unsigned>>, 10,
  //     8, 512>
  //     synt(i);
  return synt.get_value(static_cast<typename fixe_t::storage_t>(i));
}
#else


int main() {
  return 0;
}
#endif
