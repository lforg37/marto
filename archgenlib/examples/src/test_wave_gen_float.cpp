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
  static constexpr auto freq = static_cast<float>(fNum);
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
    return ul_type{}.result(time, coef);
  }
};

__attribute((always_inline)) auto test2(float i, auto& coef) {
  AdditiveSynthesizer<256, 1000> synt;
  return synt.get_value(i, coef);
}

#ifdef TARGET_VITIS
__VITIS_KERNEL auto test(float i, std::array<float, 256> coef) {
  return test2(i, coef);
}
#else

void plot() {
  std::array<float, 256> coef = {};
  coef[0] = 0.5;
  coef[1] = .5;
  for (int i = 0; i < (1 << 12); i++)
    std::cout << i << ","<< test2(std::ldexp(i, -12), coef) << std::endl;
}

int main() {
  plot();
}
#endif
