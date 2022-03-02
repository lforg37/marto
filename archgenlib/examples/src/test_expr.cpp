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

template<bitweight_t MSBWeight, bitweight_t LSBWeight, bool ISSigned> 
using fpnum_t = archgenlib::FPNumber<archgenlib::FPDim<MSBWeight, LSBWeight, ISSigned>, hint::BitIntWrapper>; 

int main() {
  int val = 17;
  archgenlib::Variable<fpnum_t<5, -4, false>> a{{val}};
  using const_valtype = hint::detail::bitint_base_t<false, 18>;
  using dim_t = archgenlib::FPDim<14, -3, false>;
  using const_t = archgenlib::FixedConstant<dim_t, const_valtype{165}>;
  archgenlib::Constant<const_t> b{};
  auto c = a + b - b;
  auto res = archgenlib::evaluate<-14>(c);
  std::cout << res << std::endl;
  return 0;
}
