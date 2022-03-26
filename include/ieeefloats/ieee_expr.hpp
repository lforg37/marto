#ifndef IEEE_EXPR_HPP
#define IEEE_EXPR_HPP

#include "ieeefloats/ieeetype.hpp"
#include "primitives/lzoc_shifter.hpp"
using hint::LZOC_shift;

#include "floatingpoint/expression.hpp"

template<unsigned int WE, unsigned int WF>
struct IEEEToFixedFormat
{
	private:
		using ieeedim = IEEEDim<WE, WF>;
		static constexpr int64_t minexp = ieeedim::MIN_NORMAL_UNBIASED_EXP - static_cast<int>(WF);
		static constexpr int64_t maxexp = ieeedim::MAX_NORMAL_UNBIASED_EXP;
	public:
		using dim = TightFixedFormat<WF, maxexp, minexp>;
};

template<unsigned int WE, unsigned int WF>
struct IEEEtoFPNum
{
	private:
		using retdim = typename IEEEToFixedFormat<WE, WF>::dim;
		template<template<unsigned int, bool> class Wrapper>
		using rettype = FixedNumber<retdim, Wrapper>;
		using ieeedim = IEEEDim<WE, WF>;

	public:
		template<template<unsigned int, bool> class Wrapper>
		static inline rettype<Wrapper> compute(IEEENumber<WE, WF, Wrapper> const & input)
		{
			auto sign = input.template get<WF+WE>();
			auto exp = input.template slice<WF+WE-1, WF>();
			auto frac = input.template slice<WF-1, 0>();
			auto expIsFullOne = exp.and_reduction();
			auto subNormal = exp.nor_reduction();
			auto fracIsFullZero = frac.nor_reduction();
			auto isNaN = expIsFullOne & fracIsFullZero.invert();
			auto isInf = expIsFullOne & fracIsFullZero;
			auto isZero = subNormal & fracIsFullZero;

			Wrapper<retdim::WE, false> bias{ieeedim::BIAS};
			Wrapper<retdim::WE, false> mbias{-ieeedim::BIAS};
			auto normalExp = exp.template leftpad<retdim::WE>().modularSub(bias).as_signed();


			auto lzocShift = LZOC_shift<WF, WF>(frac, {{0}});
			auto lzoc = lzocShift.lzoc;
			auto snFrac = lzocShift.shifted.template slice<WF - 2, 0>().template rightpad<WF>();
			auto snExp = mbias.modularSub(lzoc.template leftpad<retdim::WE>()).as_signed();

			auto rExp = Wrapper<retdim::WE, true>::mux(subNormal, snExp, normalExp);
			auto rFrac = Wrapper<WF, false>::mux(subNormal, snFrac, frac);

			return {rFrac, rExp, sign, isInf, isNaN, isZero};
		}
};

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
FPExpr<typename IEEEToFixedFormat<WE, WF>::type, Wrapper> inline to_expr(IEEENumber<WE, WF, Wrapper> const in)
{
	return to_expr(IEEEtoFPNum<WE, WF>::compute(in));
}
#endif // IEEE_EXPR_HPP
