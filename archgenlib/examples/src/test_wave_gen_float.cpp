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
  static constexpr float wfreq = fNum * 2. * std::numbers::pi_v<double> ;
  wave_gen(float phase) : phase{phase} {}

  auto get_from_time() {
    auto res = std::sin(phase);
    phase += wfreq;
    return res;
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
  osc o{0.};

  auto result(std::array<float, NBOsc>& coef) {
    return coef[start_idx - 1] * o.get_from_time();
  }
};

template<unsigned start_idx, unsigned last_idx, unsigned int NBOsc, auto max_frequency>
struct OscillatorBench<start_idx, last_idx, NBOsc, max_frequency, std::enable_if_t<start_idx < last_idx>>
{
  static constexpr unsigned mid_point = (last_idx + start_idx) / 2;
  using low_type = OscillatorBench<start_idx, mid_point, NBOsc, max_frequency>;
  using high_type = OscillatorBench<mid_point + 1, last_idx, NBOsc, max_frequency>;
  low_type l;
  high_type h;
  auto result(auto& coef) {
    auto res_low = l.result(coef);
    return res_low + h.result(coef);
  }
};

}

template<unsigned int NBOsc, auto max_frequency>
struct AdditiveSynthesizer {
  using ul_type = detail::OscillatorBench<1, NBOsc, NBOsc, max_frequency>;
  ul_type u;
  auto get_value(auto& coef) {
    return u.result(coef);
  }
};

template<unsigned int NBOsc, auto max_frequency>
AdditiveSynthesizer<NBOsc, max_frequency> additive_synth{};

__attribute((always_inline)) auto test2(auto& coef) {
  auto& synt = additive_synth<256, 1000>;
  return synt.get_value(coef);
}

#ifdef __VITIS_KERNEL
__VITIS_KERNEL auto test(std::array<float, 256> coef) {
  return test2(coef);
}
#else

void plot() {
  std::array<float, 256> coef = {};
  coef[0] = 0.5;
  coef[1] = .5;
  for (int i = 0; i < (1 << 12); i++)
    std::cout << i << ","<< test2(coef) << std::endl;
}

int main() {
  plot();
}
#endif
