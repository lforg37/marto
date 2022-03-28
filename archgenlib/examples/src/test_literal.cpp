#include "fixedpoint/literal.hpp"
#include <type_traits>

int main() {
  static_assert(
      std::is_same<decltype(0x0_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<0, 0, unsigned>, 0>>>::value,
      "");
  static_assert(
      std::is_same<decltype(0x1_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<1, 0, unsigned>, 1>>>::value,
      "");
  static_assert(
      std::is_same<decltype(0x2_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<2, 0, unsigned>, 2>>>::value,
      "");
  static_assert(
      std::is_same<decltype(0x3_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<2, 0, unsigned>, 3>>>::value,
      "");
  static_assert(
      std::is_same<decltype(0x1.1p0_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<1, -4, unsigned>, 17>>>::value,
      "");
  static_assert(
      std::is_same<decltype(0x1.8p0_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<1, -1, unsigned>, 3>>>::value,
      "");
  static_assert(
      std::is_same<decltype(0x8.8p0_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<4, -1, unsigned>, 17>>>::value,
      "");
  static_assert(
      std::is_same<decltype(0x008.800p0_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<4, -1, unsigned>, 17>>>::value,
      "");
  static_assert(
      std::is_same<decltype(0x1.8p-2_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<-1, -3, unsigned>, 3>>>::value,
      "");
  static_assert(
      std::is_same<decltype(0x1000ff0f0f0f0f0f0f0f0.0p0_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<81, 0, unsigned>,
                       ((unsigned _BitInt(81))0x1000f << 64) |
                           (unsigned _BitInt(81))0xf0f0f0f0f0f0f0f0>>>::value,
      "");
  static_assert(
      std::is_same<decltype(0x1000ff0f0f0f0f0f0f0f0.0p-10_cst),
                   archgenlib::Constant<archgenlib::FixedConstant<
                       archgenlib::FixedFormat<71, -10, unsigned>,
                       ((unsigned _BitInt(81))0x1000f << 64) |
                           (unsigned _BitInt(81))0xf0f0f0f0f0f0f0f0>>>::value,
      "");
}
