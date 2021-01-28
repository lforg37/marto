#pragma once
#include <cstdio>
//#include <tools/printing.hpp>

#include "posit_dim.hpp"
#include "posit_in_place_round.hpp"

#ifdef POSIT_VALPROD_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif


template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositProd<N, WES, Wrapper> PositIF_to_PositProd(PositIntermediateFormat<N, WES, Wrapper, true> val)
{
	auto val_frac = val.getSignificand();
	auto significand = val_frac.concatenate(
			Wrapper<PositDim<N, WES>::WF+1, false>::generateSequence({0})
		);


	auto expt = val.getExp().as_signed().template leftpad<PositDim<N, WES>::ProdExpSize>().as_unsigned();

	return {val.getIsNaR(), expt, val.getSignBit(), significand};
}

//TODO adapt exact
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
	auto exp = val.getExp();

	// hint<1> isZero = not(((hint<4>) val.getSignificand().slice(PositDim<N>::ProdSignificandSize - 1, PositDim<N>::ProdSignificandSize - 4)).or_reduce());
	auto isZero = sign.invert() & implicitBit.invert();

	//cerr << to_string(isZero) << endl;
	//cerr << to_string(exp) << endl;

	// exp.print();
	auto emax = Wrapper<PROD_WE, false>{PositDim<N, WES>::EMax};
	auto exp_pos = exp.template get<PROD_WE - 1>().invert(); // Exponent sign bit
	auto eabs_mask = Wrapper<PROD_WE, false>::generateSequence(exp_pos);
	auto exp_aabs = exp ^ eabs_mask; // exp_aabs is either -|exp| or -|exp + 1|

	// Get sign bit of Emax - |abs| to check if exponent overflow
	auto expOverFlow = emax.addWithCarry(exp_aabs, exp_pos).template get<PROD_WE - 1>() & isNaR.invert() & isZero.invert();


	auto isMinPos = expOverFlow & exp_pos.invert();
	auto isMaxPos = expOverFlow & exp_pos;

	//cerr << to_string(isMaxPos) << endl;
	//cerr << to_string(isMinPos) << endl;

	auto minposval = Wrapper<PositDim<N, WES>::PIFSize, false>::mux(
					sign,
					PositIntermediateFormat<N, WES, Wrapper, true>::getMinNeg(),
					PositIntermediateFormat<N, WES, Wrapper, true>::getMinPos()
				);

	auto maxposval = Wrapper<PositDim<N, WES>::PIFSize, false>::mux(
				sign,
				PositIntermediateFormat<N, WES, Wrapper, true>::getMaxNeg(),
				PositIntermediateFormat<N, WES, Wrapper, true>::getMaxPos()
			);

	auto specialval = Wrapper<PositDim<N, WES>::PIFSize, false>::mux(
				isMaxPos,
				maxposval,
				minposval
			);

	auto ret = Wrapper<PositDim<N, WES>::UPIFSize, false>::mux(
				expOverFlow,
				specialval.template leftpad<PositDim<N, WES>::UPIFSize>(),
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

#ifdef POSIT_VALPROD_DEBUG
	cerr << "=== POSIT_PROD_TO_POSIT_IF" << endl;
	cerr << "isNaR: " << to_string(isNaR) << endl;
	cerr << "fraction: " << to_string(fraction) << endl;
	cerr << "implicitBit: " << to_string(implicitBit) << endl;
	cerr << "sign: " << to_string(sign) << endl;
	cerr << "resultFraction: " << to_string(resultFraction) << endl;
	cerr << "resultGuardBit: " << to_string(resultGuardBit) << endl;
	cerr << "resultStickyBit: " << to_string(resultStickyBit) << endl;
	cerr << "exp: " << to_string(exp) << endl;
	cerr << "isZero: " << to_string(isZero) << endl;
	cerr << "emax: " << to_string(emax) << endl;
	cerr << "exp_pos: " << to_string(exp_pos) << endl;
	cerr << "eabs_mask: " << to_string(eabs_mask) << endl;
	cerr << "exp_aabs: " << to_string(exp_aabs) << endl;
	cerr << "expOverFlow: " << to_string(expOverFlow) << endl;
	cerr << "isMinPos: " << to_string(isMinPos) << endl;
	cerr << "isMaxPos: " << to_string(isMaxPos) << endl;
	cerr << "minposval: " << to_string(minposval) << endl;
	cerr << "maxposval: " << to_string(maxposval) << endl;
	cerr << "specialval: " << to_string(specialval) << endl;
	cerr << "ret: " << to_string(ret) << endl;
	cerr << "=========================" << endl;
#endif
	return ret;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, true> PositProd_to_PositIF_in_place_rounding(PositProd<N, WES, Wrapper> const & val)
{
	auto prod_inexact = PositProd_to_PositIF(val);
	auto ret = in_place_rounder(prod_inexact);
	return ret;
}
