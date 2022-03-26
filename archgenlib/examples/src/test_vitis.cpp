#include <cstdint>
#include <string_view>
#include <type_traits>

#include "bitint_tools/type_helpers.hpp"
#include "hint.hpp"

#include "fixedpoint/expression.hpp"
#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"
#include "fixedpoint/literal.hpp"

using archgenlib::bitweight_t;  

constexpr bitweight_t outprec = -10;

using fpdim_t = archgenlib::FixedFormat<5, outprec, unsigned>;
using fpnum_t = archgenlib::FixedNumber<fpdim_t>;

using storage_t =
    hint::detail::bitint_base_t<fpdim_t::is_signed, fpdim_t::width>;

#ifndef __VITIS_KERNEL
#warning "unsupported compiler"
#define __VITIS_KERNEL
#endif

auto test(int i) {
  auto val = static_cast<storage_t>(i);
  archgenlib::Variable<fpnum_t> a{{val}};
  using const_valtype = hint::detail::bitint_base_t<false, 16>;
  using dim_t = archgenlib::FixedFormat<14, -1, unsigned>;
  using const_t = archgenlib::FixedConstant<dim_t, const_valtype{3}>;
  archgenlib::Constant<const_t> b{};
  auto c = archgenlib::sin(a) * b;
  return archgenlib::evaluate<outprec>(c);
}

__VITIS_KERNEL auto
function_sin(archgenlib::FixedNumber<archgenlib::FixedFormat<5, -12, unsigned>> val) {
  archgenlib::Variable x{val};
  auto c = 0x3.p-1_cst;
  auto expr = sin(x + c);
  return archgenlib::evaluate<-17>(expr);
}
