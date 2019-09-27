#ifndef IEEE_MULTIPLIER_HPP
#define IEEE_MULTIPLIER_HPP

#include "primitives/lzoc.hpp"
#include "tools/static_math.hpp"

#include "ieeetype.hpp"

namespace hint {
	signed int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
	IEEENumber<WE, WF, Wrapper> ieee_product(IEEENumber<WE, WF, Wrapper> i0, IEEENumber<WE, WF, Wrapper> i1)
	{
		auto exp0 = i0.getExponent();
		auto exp1 = i1.getExponent();

		auto frac0 = i0.getFractionnalPart();
		auto frac1 = i1.getFractionnalPart();

		auto unsorted = exp0 < exp1;

		auto expMax = Wrapper<WE, false>::mux(unsorted, exp1, exp0);
		auto expMin = Wrapper<WE, false>::mux(unsorted, exp0, exp1);

		auto fracMax = Wrapper<WF, false>::mux(unsorted, frac1, frac0);
		auto fracMin = Wrapper<WF, false>::mux(unsorted, frac0, frac1);

		auto lzocMin = hint::lzoc_wrapper(fracMin);

		//---------- Special cases detection -----------------------
		auto isSubnormalMax = expMax.or_reduction().invert();
		auto isSubnormalMin = expMin.or_reduction().invert();

		auto isMaxExpMax = expMax.and_reduction();
		auto isMaxExpMin = expMin.and_reduction();

		auto isFracNotEmptyMax = fracMax.or_reduction();
		auto isFracNotEmptyMin = fracMin.or_reduction();

		auto isInfiniteMax = isMaxExpMax & (isFracNotEmptyMax.invert());
		auto isInfiniteMin = isMaxExpMin & (isFracNotEmptyMin.invert());

		auto isNaNMax = isMaxExpMax & isFracNotEmptyMax;
		auto isNaNMin = isMaxExpMin & isFracNotEmptyMin;

		auto isEqToZeroMax = isSubnormalMax & (isFracNotEmptyMax.invert());
		auto isEqToZeroMin = isSubnormalMin & (isFracNotEmptyMin.invert());

		auto oneIsNaN = isNaNMax | isNaNMin;
		auto oneIsInfty = isInfiniteMax | isInfiniteMin;
		auto oneIsZero = isEqToZeroMax | isEqToZeroMin;

		auto isResNaN = oneIsNaN | (oneIsInfty & oneIsZero);
		auto isResZero = oneIsZero & (oneIsNaN.invert()) & (oneIsInfty.invert());

		//---------------------------------------------------------
		auto expOffset0 = Wrapper<hint::Static_Val<WE>::_storage, false>::mux(
					isSubnormalMax,
					lzocMaxz,
					{0}
				);
		auto expOffset1 = Wrapper<hint::Static_Val<WE>::_storage, false>::mux(
					isSubnormal0,
					lzoc0,
					{0}
				);

		auto o_sign = i0.getSign() ^ i1.getSign();
		//-------------------- Handling exponents ----------------------

		auto exp = exp0 + exp1;
		auto prod = (isSubnormal0.invert().concat(frac0)) * (isSubnormal1.invert().concat(frac1));

		constexpr size_t prodsize = hint::Arithmetic_Prop<WF+1, WF+1>::_prodSize;

		auto high_frac = prod.template slice<prodsize - 2, prodsize - (WE + 1)>();
		auto low_frac = prod.template slice<prodsize - 3, prodsize - (WE + 2)>();


	}
}
template<un
#endif // IEEE_MULTIPLIER_HPP
