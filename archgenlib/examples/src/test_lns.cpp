#include <cstdint>
#include <string_view>
#include <type_traits>

#include "bitint_tools/type_helpers.hpp"
#include "fixedpoint/evaluator.hpp"
#include "fixedpoint/functions.hpp"
#include "hint.hpp"

#include "fixedpoint/expression.hpp"
#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"

#ifndef __VITIS_KERNEL
#warning "unsupported compiler"
#define __VITIS_KERNEL
#endif

template<typename storage_t>
class LNSNumber;

template <typename> static constexpr bool is_lns_number_v = false;
template <typename T>
static constexpr bool is_lns_number_v<LNSNumber<T>> = true;

template<typename T>
concept is_lns_number = is_lns_number_v<T>;

// nan, sign, FPBNumber
template<typename storage_t>
class LNSNumber {
  using dim_t = typename storage_t::dimension_t;
  template <is_lns_number Other>
  using add_bitwidth_plus1 = LNSNumber<archgenlib::FixedNumber<archgenlib::FixedFormat<
      std::max(dim_t::msb_weight, Other::dim_t::msb_weight) + 1,
      std::min(dim_t::lsb_weight, Other::dim_t::lsb_weight), unsigned>>>;

  storage_t value_;
  bool is_nan;
  bool is_zero;
  bool is_neg;
  public:
  static LNSNumber nan() {
    return {0, true, false, false};
  }
  static LNSNumber zero() {
    return {0, true, false, false};
  }
  constexpr LNSNumber(storage_t const &val) : value_{val} {}
  storage_t value() const { return value_; }
  auto as_hint() const { return value_.as_hint(); }

  template <is_lns_number OtherTy>
  add_bitwidth_plus1<OtherTy> operator*(const OtherTy &other) {
    return {value() + other.value()};
  }
  template<is_lns_number OtherTy>
  add_bitwidth_plus1<OtherTy> operator/(const OtherTy& other) {
    return {value() - other.value()};
  }
  template <is_lns_number OtherTy>
  add_bitwidth_plus1<OtherTy> operator+(const OtherTy &other){
    archgenlib::Variable<add_bitwidth_plus1<OtherTy>> var{};
    return {archgenlib::evaluate<0>(archgenlib::log2(var)).value()};
  }
};

using archgenlib::bitweight_t;  

constexpr bitweight_t outprec = -10;

using fpdim_t = archgenlib::FixedFormat<5, outprec, unsigned>;
using fpnum_t = archgenlib::FixedNumber<fpdim_t>;

using storage_t =
    hint::detail::bitint_base_t<fpdim_t::is_signed, fpdim_t::width>;

__VITIS_KERNEL auto
function_sin(LNSNumber<fpnum_t> val, LNSNumber<fpnum_t> val2, LNSNumber<fpnum_t> val3) {
  return val + val2 * val3;
}
