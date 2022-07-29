#include <compare>
#include <cstdint>
#include <string_view>
#include <type_traits>

#include "bitint_tools/type_helpers.hpp"
#include "fixedpoint/header_gen_evaluator.hpp"
#include "hint.hpp"

#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"
#include "fixedpoint/literal.hpp"

#ifndef __VITIS_KERNEL
#warning "unsupported compiler"
#define __VITIS_KERNEL
#endif

template<int base, typename storage_t>
class LNSNumber;

template <typename> static constexpr bool is_lns_number_v = false;
template <int base, typename T>
static constexpr bool is_lns_number_v<LNSNumber<base, T>> = true;

template<typename T>
concept is_lns_number = is_lns_number_v<T>;

template <typename T, typename T2>
concept is_compatible_lns =
    is_lns_number_v<T> &&is_lns_number_v<T2> &&T::base == T2::base;

template<int ba, typename storage_t>
class LNSNumber {
  template<int b, typename s>
  friend class LNSNumber;
public:
  static constexpr int base = ba;
  using dim_t = typename storage_t::format_t;
private:
  template <is_lns_number Other>
  using add_bitwidth_plus1 = archgenlib::FixedNumber<archgenlib::FixedFormat<
      std::max(dim_t::msb_weight, Other::dim_t::msb_weight) + 1,
      std::min(dim_t::lsb_weight, Other::dim_t::lsb_weight), unsigned>>;
  using base_cst =
      archgenlib::Constant<
      archgenlib::FixedConstant<archgenlib::FixedFormat<10, 0, signed>, base>>;

  storage_t value_ = 0;
  bool is_nan_ = 0;
  bool is_neg_ = 0;
  void normalize() {
    if (value_.value() == 0 && is_neg_)
      is_neg_ = false;
  }

  LNSNumber(storage_t s, bool nan, bool neg) : value_{s}, is_nan_{nan}, is_neg_{neg} {}

  public:
  static LNSNumber nan() {
    return {0, true, false};
  }
  static LNSNumber zero() {
    return {0, false, false};
  }

  bool is_nan() {
    return is_nan_;
  }
  bool is_neg() {
    return is_neg_;
  }
  bool is_zero() {
    return !is_nan() && value_.value() == 0;
  }
  storage_t value() { return value_; }

  constexpr LNSNumber(storage_t const &val) : value_{val} {}
  storage_t value() const { return value_; }
  auto as_hint() const { return value_.as_hint(); }

  template <is_compatible_lns<LNSNumber> OtherTy>
  LNSNumber<base, add_bitwidth_plus1<OtherTy>> operator*(const OtherTy &other) {
    LNSNumber<base, add_bitwidth_plus1<OtherTy>> res{value() + other.value()};
    res.is_nan_ = is_nan_ || other.is_nan_;
    res.is_neg_ = is_neg_ ^ other.is_neg_;
    res.normalize();
    return res;
  }
  template<is_compatible_lns<LNSNumber> OtherTy>
  LNSNumber<base, add_bitwidth_plus1<OtherTy>> operator/(const OtherTy& other) {
    LNSNumber<base, add_bitwidth_plus1<OtherTy>> res{value() - other.value()};
    res.is_nan_ = is_nan_ || other.is_nan_ || other.is_zero();
    res.is_neg_ = is_neg_ ^ other.is_neg_;
    res.normalize();
    return res;
  }
  template <is_compatible_lns<LNSNumber> OtherTy>
  std::partial_ordering operator<=>(OtherTy Other) {
    if (Other.is_nan_ || is_nan_)
      return std::partial_ordering::unordered;
    if (Other.is_neg_ != is_neg_) {
      if (is_neg_)
        return std::partial_ordering::less;
      return std::partial_ordering::greater;
    }
    if (is_neg_)
      return Other.value_.operator<=>(value_);
    return value_.operator<=>(Other.value_);
  }
  template <is_compatible_lns<LNSNumber> OtherTy>
  LNSNumber<base, add_bitwidth_plus1<OtherTy>> operator+(const OtherTy &other) {
    archgenlib::Variable<add_bitwidth_plus1<OtherTy>> var{other.value() -
                                                          value()};
    auto res = LNSNumber<base, add_bitwidth_plus1<OtherTy>>::zero();
    res.value_ =
        (value() +
         archgenlib::evaluate<0>(
             archgenlib::log(0x1.p0_cst + archgenlib::pow(base_cst{}, var)) /
             archgenlib::log(base_cst{})))
            .value();
    res.is_neg_ =
        (value().value() >= other.value().value()) ? is_neg_ : other.is_neg_;
    res.is_nan_ = is_nan_ || other.is_nan_;
    res.normalize();
    return res;
  }
  template <is_compatible_lns<LNSNumber> OtherTy>
  LNSNumber<base, add_bitwidth_plus1<OtherTy>> operator-(const OtherTy &other) {
    archgenlib::Variable<add_bitwidth_plus1<OtherTy>> var{other.value() -
                                                          value()};
    auto res = LNSNumber<base, add_bitwidth_plus1<OtherTy>>::zero();
    res.value_ =
        (value() + archgenlib::evaluate<0>(
                       archgenlib::log(archgenlib::abs(
                           0x1.p0_cst - archgenlib::pow(base_cst{}, var))) /
                       archgenlib::log(base_cst{})))
            .value();
    res.is_neg_ =
        (value().value() >= other.value().value()) ? is_neg_ : !other.is_neg_;
    res.is_nan_ = is_nan_ || other.is_nan_;
    res.normalize();
    return res;
  }
  LNSNumber operator-() {
    is_neg_ = !is_neg_;
    normalize();
  }
};

using archgenlib::bitweight_t;  

constexpr bitweight_t outprec = -10;

using fpdim_t = archgenlib::FixedFormat<5, outprec, unsigned>;
using fpnum_t = archgenlib::FixedNumber<fpdim_t>;

using storage_t =
    hint::detail::bitint_base_t<fpdim_t::is_signed, fpdim_t::width>;

__VITIS_KERNEL auto
function_sin(LNSNumber<2, fpnum_t> val, LNSNumber<2, fpnum_t> val2, LNSNumber<2, fpnum_t> val3, LNSNumber<2, fpnum_t> val4) {
  return val + val2 * val3 - val4;
}
