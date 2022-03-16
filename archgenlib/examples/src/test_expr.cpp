#include <iostream>
#include <type_traits>

#include "hint.hpp"

#include "fixedpoint/expression.hpp"
#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"

#ifdef INCLUDE_GENERATED_HEADER
#include "specialization_header.hpp"
#endif

template <typename T> static constexpr bool ok = false;

template <archgenlib::ExpressionType T> static constexpr bool ok<T> = true;

using archgenlib::bitweight_t;

template <bitweight_t MSBWeight, bitweight_t LSBWeight, bool ISSigned>
using fpnum_t =
    archgenlib::FPNumber<archgenlib::FPDim<MSBWeight, LSBWeight, ISSigned>>;

int main() {
  unsigned _BitInt(10) val = 17;
  archgenlib::Variable<fpnum_t<5, -4, false>> a{{val}};
  using const_valtype = hint::detail::bitint_base_t<false, 16>;
  using dim_t = archgenlib::FPDim<14, -1, false>;
  using const_t = archgenlib::FixedConstant<dim_t, const_valtype{3}>;
  archgenlib::Constant<const_t> b{};
  auto c = sin(a) * b;
  auto res = archgenlib::evaluate<-4>(c);
  std::cout << static_cast<int>(res.value()) << std::endl;
  return 0;
}
