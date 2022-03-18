#ifndef FIXEDPOINT_FIXEDPOINT_HPP
#define FIXEDPOINT_FIXEDPOINT_HPP

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <limits>
#include <string_view>
#include <type_traits>

#ifndef BITINT_BACKEND
#define BITINT_BACKEND
#endif
#include "hint.hpp"
#include "bitint_tools/bitint_constant.hpp"
#include "bitint_tools/type_helpers.hpp"

namespace archgenlib {
using bitweight_t = std::int32_t;
using vecwidth_t = std::uint32_t;

template <bitweight_t MSBWeight, bitweight_t LSBWeight, bool IsSigned>
struct FPDim {
  static_assert(MSBWeight >= LSBWeight,
                "MSB weight cannot be smaller than LSBs");

  /**
   * @brief msb_weight Weight of the most significant bit
   */
  static constexpr bitweight_t msb_weight = MSBWeight;

  /**
   * @brief lsb_weight
   */
  static constexpr bitweight_t lsb_weight = LSBWeight;

  /**
   * @brief WFF Full significand width
   */
  static constexpr vecwidth_t width = MSBWeight - LSBWeight + 1;

  /**
   * @brief is_signed is the number signed
   */
  static constexpr bool is_signed = IsSigned;
};

namespace detail {
template <typename T> constexpr bool is_fpdim = false;

template <bitweight_t MSBWeight, bitweight_t LSBWeight, bool IsSigned>
constexpr bool is_fpdim<FPDim<MSBWeight, LSBWeight, IsSigned>> = true;
} // namespace detail

template <typename T>
concept FPDimType = detail::is_fpdim<T>;

template <std::integral IT>
using fpdim_from_integral =
    FPDim<std::numeric_limits<IT>::digits + std::is_signed_v<IT>, 0,
          std::is_signed_v<IT>>;

template <FPDimType Dim> class FPNumber {
  using storage_t = hint::detail::bitint_base_t<Dim::is_signed, Dim::width>;

  storage_t value_;

public:
  using dimension_t = Dim;
  static constexpr auto width = Dim::width;
  constexpr FPNumber(storage_t const &val) : value_{val} {}
  storage_t value() const { return value_; }
  hint::BitIntWrapper<Dim::width, Dim::is_signed> as_hint() const {
    return {value_};
  }
};

namespace detail {
template <typename T> constexpr bool is_fp_num = false;

template <FPDimType Dim> constexpr bool is_fp_num<FPNumber<Dim>> = true;
} // namespace detail

template <typename T>
concept FPNumberType = detail::is_fp_num<T>;

// FixedPointConstant
template <typename T, typename Dim>
concept FixedBitIntCst = hint::BitIntConstT<T> && FPDimType<Dim> && T::width ==
Dim::width &&T::isSigned == Dim::is_signed;

template <FPDimType Dim,
          hint::detail::bitint_base_t<Dim::is_signed, Dim::width> Val>
struct FixedConstant {
  using dimension_t = Dim;
  using value_type = hint::detail::bitint_base_t<Dim::is_signed, Dim::width>;
  static constexpr value_type value = Val;
};

namespace detail {
template <typename T> constexpr bool _is_fixed_constant = false;
template <FPDimType Dim,
          hint::detail::bitint_base_t<Dim::is_signed, Dim::width> Val>
constexpr bool _is_fixed_constant<FixedConstant<Dim, Val>> = true;
} // namespace detail

template <typename T>
concept FixedConstantType = detail::_is_fixed_constant<T>;

namespace detail {}

} // namespace archgenlib

#endif