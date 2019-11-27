#pragma once

#include <cstdio>
#include <tools/printing.hpp>

#include "posit_dim.hpp"

using hint::to_string;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositProd<N, WES, Wrapper> PositIF_to_PositProd(PositIntermediateFormat<N, WES, Wrapper, true> val)
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
inline PositIntermediateFormat<N, WES, Wrapper, false> PositProd_to_PositIF(PositProd<N, WES, Wrapper> val)
{
	auto isNaR = val.getIsNaR();
	auto fraction = val.getSignificand();
	auto implicitBit = fraction.template get<PositDim<N, WES>::ProdSignificandSize-1>();
	auto sign = val.getSignBit();
	constexpr unsigned int WF = PositDim<N, WES>::WF;
	constexpr unsigned int PROD_WF = PositDim<N, WES>::ProdSignificandSize;
	constexpr unsigned int PROD_WE = PositDim<N, WES>::ProdExpSize;

	auto resultFraction = fraction.template slice<PROD_WF - 2, PROD_WF-1-WF>();
	//cerr << to_string(resultFraction) << endl;

	auto resultGuardBit = fraction.template get<PROD_WF-2-WF>();

	auto resultStickyBit = fraction.template slice<PROD_WF-3 -WF, 0>().or_reduction();
	auto expt = val.getExp();

	// hint<1> isZero = not(((hint<4>) val.getSignificand().slice(PositDim<N>::ProdSignificandSize - 1, PositDim<N>::ProdSignificandSize - 4)).or_reduce());
	auto isZero = sign.invert() & implicitBit.invert() & expt.or_reduction().invert();

	auto exp = Wrapper<PROD_WE + 1, false>::mux(
					isZero,
					{0},
					expt.template leftpad<PROD_WE + 1>().modularSub({PositDim<N, WES>::EXP_BIAS-1})
				);
	//cerr << to_string(isZero) << endl;
	//cerr << to_string(exp) << endl;

	// exp.print();
	auto maxPosCheck = exp.modularSub({2*PositDim<N, WES>::EXP_BIAS})
				.template get<PROD_WE>();
	auto isMinPos = exp.template get<PositDim<N, WES>::ProdExpSize>() & isNaR.invert();
	auto isMaxPos = maxPosCheck.invert() & isNaR.invert();

	//cerr << to_string(isMaxPos) << endl;
	//cerr << to_string(isMinPos) << endl;

	auto minposval = Wrapper<PositDim<N, WES>::ValSize, false>::mux(
					sign,
					PositIntermediateFormat<N, WES, Wrapper, true>::getMinNeg(),
					PositIntermediateFormat<N, WES, Wrapper, true>::getMinPos()
				);

	auto maxposval = Wrapper<PositDim<N, WES>::ValSize, false>::mux(
				sign,
				PositIntermediateFormat<N, WES, Wrapper, true>::getMaxNeg(),
				PositIntermediateFormat<N, WES, Wrapper, true>::getMaxPos()
			);

	auto specialval = Wrapper<PositDim<N, WES>::ValSize, false>::mux(
				isMaxPos,
				maxposval,
				minposval
			);

	auto ret = Wrapper<PositDim<N, WES>::ValSize, false>::mux(
				isMaxPos.bitwise_or(isMinPos),
				specialval,
				PositIntermediateFormat<N, WES, Wrapper, false>(
								resultGuardBit,
								resultStickyBit,
								isNaR,
								exp.template slice<PositDim<N, WES>::WE - 1, 0>(), //Warning : biased exp
								val.getSignBit(),
								implicitBit,
								resultFraction
							)
			);
	/*
	cerr << to_string(resultGuardBit) << endl;
	cerr << to_string(resultStickyBit) << endl;
	cerr << to_string(isNaR) << endl;
	cerr << to_string(exp) << endl;
	cerr << to_string(val.getSignBit()) << endl;
	cerr << to_string(implicitBit) << endl;
	cerr << to_string(resultFraction) << endl;
	cerr << to_string(ret) << endl;
	*/
	return ret;
}
