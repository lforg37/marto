#pragma once

#include <cstdio>

#include "posit_dim.hpp"

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositProd<N, WES, Wrapper> PositIF_to_PositProd(PositIntermediateFormat<N, WES, Wrapper> val)
{
	auto val_frac = val.getSignificand();
	auto significand = val_frac.concatenate(
			Wrapper<PositDim<N, WES>::WF+1, false>::generateSequence({0})
		);


	auto expt = val.getExp().template leftpad<PositDim<N, WES>::ProdExpSize>();

	auto exponent = Wrapper<PositDim<N, WES>::ProdExpSize, false>::mux(
					val.isZero(),
					{0},
					expt.modularAdd(Wrapper<PositDim<N, WES>::ProdExpSize, false>{PositDim<N, WES>::EXP_BIAS - 1})
				);

	return {val.getIsNaR(), exponent, val.getSignBit(), significand};
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper> PositProd_to_PositIF(PositProd<N, WES, Wrapper> val)
{
	auto isNaR = val.getIsNaR();
	auto fraction = val.getSignificand();
	auto implicitBit = fraction.template get<PositDim<N, WES>::ProdSignificandSize-1>();
	auto sign = val.getSignBit();
	constexpr unsigned int WF = PositDim<N, WES>::WF;
	constexpr unsigned int prod_WF = PositDim<N, WES>::ProdSignificandSize;

	auto resultFraction = fraction.template slice<prod_WF - 2, prod_WF-1-WF>();

	auto resultGuardBit = fraction.template get<prod_WF-2-WF>();

	auto resultStickyBit = fraction.template slice<prod_WF-3 -WF, 0>().or_reduction().invert();
	auto expt = val.getExp();

	// hint<1> isZero = not(((hint<4>) val.getSignificand().slice(PositDim<N>::ProdSignificandSize - 1, PositDim<N>::ProdSignificandSize - 4)).or_reduce());
	auto isZero = val.getSignBit().invert().bitwise_and(implicitBit.invert()).bitwise_and(expt.or_reduction().invert());

	auto exp = Wrapper<PositDim<N, WES>::ProdExpSize + 1, false>::mux(
					isZero,
					{0},
					expt.template leftpad<PositDim<N, WES>::ProdExpSize + 1>().modularSub({PositDim<N, WES>::EXP_BIAS-1})
				);

	// exp.print();
	auto maxPosCheck = exp.modularSub({2*PositDim<N, WES>::EXP_BIAS});
	auto isMinPos = exp.template get<PositDim<N, WES>::ProdExpSize>();
	auto isMaxPos = maxPosCheck.template get<PositDim<N, WES>::ProdExpSize>().invert().bitwise_and(isNaR.invert());

	// fprintf(stderr, "isminpos: ");
	// isMinPos.print();
	// fprintf(stderr, "ismaxpos: ");
	// isMaxPos.print();

	auto minposval = Wrapper<PositDim<N, WES>::ValSize, false>::mux(
					sign,
					PositIntermediateFormat<N, WES, Wrapper>::getMinNeg(),
					PositIntermediateFormat<N, WES, Wrapper>::getMinPos()
				);

	auto maxposval = Wrapper<PositDim<N, WES>::ValSize, false>::mux(
				sign,
				PositIntermediateFormat<N, WES, Wrapper>::getMaxNeg(),
				PositIntermediateFormat<N, WES, Wrapper>::getMaxPos()
			);

	auto specialval = Wrapper<PositDim<N, WES>::ValSize, false>::mux(
				isMaxPos,
				maxposval,
				minposval
			);

	return Wrapper<PositDim<N, WES>::ValSize, false>::mux(
				isMaxPos.bitwise_or(isMinPos),
				specialval,
				PositIntermediateFormat<N, WES, Wrapper>(
								resultGuardBit,
								resultStickyBit,
								isNaR,
								exp.template slice<PositDim<N, WES>::WE - 1, 0>(), //Warning : biased exp
								val.getSignBit(),
								implicitBit,
								resultFraction
							)
			);
}
