#ifndef FIXEDPOINT_FIXEDPOINT_HPP
#define FIXEDPOINT_FIXEDPOINT_HPP


#include <concepts>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "bitint_tools/bitint_constant.hpp"

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

template<std::integral IT>
using fpdim_from_integral = FPDim<std::numeric_limits<IT>::digits + std::is_signed_v<IT>, 0, std::is_signed_v<IT>>;

template <FPDimType Dim, template <unsigned int, bool> class Wrapper>
class FPNumber {
  using storage_t = Wrapper<Dim::width, Dim::is_signed>;

  storage_t value_;

public:
  using dimension_t = Dim;
  FPNumber(storage_t const &val) : value_{val} {}
  storage_t value() const { return value_;}
};

namespace detail {
template <typename T> constexpr bool is_fp_num = false;

template <FPDimType Dim, template <unsigned int, bool> class Wrapper>
constexpr bool is_fp_num<FPNumber<Dim, Wrapper>> = true;
} // namespace detail

template <typename T>
concept FPNumberType = detail::is_fp_num<T>;

// FixedPointConstant
template<typename T, typename Dim>
concept FixedBitIntCst = hint::BitIntConstT<T> && FPDimType<Dim> && T::width == Dim::width && T::isSigned == Dim::is_signed;

template<FPDimType Dim, hint::detail::bitint_base_t<Dim::is_signed, Dim::width> Val>
struct FixedConstant {
  using dimension_t = Dim;
  using value_type = hint::detail::bitint_base_t<Dim::is_signed, Dim::width>;
  static constexpr value_type value = Val; 
};

namespace detail {
  template<typename T> constexpr bool _is_fixed_constant = false;
  template<FPDimType Dim, hint::detail::bitint_base_t<Dim::is_signed, Dim::width> Val> 
  constexpr bool _is_fixed_constant<FixedConstant<Dim, Val>> = true;
}

template<typename T>
concept FixedConstantType = detail::_is_fixed_constant<T>;

} // namespace archgenlib

#endif
