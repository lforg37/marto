#ifndef FIXEDPOINT_FIXEDPOINT_HPP
#define FIXEDPOINT_FIXEDPOINT_HPP


#include <concepts>
#include <cstdint>
#include <limits>
#include <type_traits>

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
  using dim_t = Dim;
  using storage_t = Wrapper<Dim::width, Dim::is_signed>;

  storage_t value_;

public:
  FPNumber(storage_t const &val) : value_{val} {}
};

namespace detail {
template <typename T> constexpr bool is_fp_num = false;

template <FPDimType Dim, template <unsigned int, bool> class Wrapper>
constexpr bool is_fp_num<FPNumber<Dim, Wrapper>> = true;
} // namespace detail

template <typename T>
concept FPNumberType = detail::is_fp_num<T>;

} // namespace archgenlib

#endif
