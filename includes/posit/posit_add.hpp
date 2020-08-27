#pragma once
#include <cstdio>
#include "tools/static_math.hpp"
#include "primitives/shifter_sticky.hpp"
#include "primitives/lzoc_shifter.hpp"
#include "primitives/table.hpp"
#include "helpers/bit_sequence_generator.hpp"

#include "posit_dim.hpp"
#include "posit_in_place_round.hpp"
#ifdef POSIT_ADDER_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif

/*
	Uppon testing, the cheapest way to perform the add_sub component is to
	negate the input before the addition instead of merging the negation in
	the operator.
*/

using namespace hint;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, false> posit_add(
		PositIntermediateFormat<N, WES, Wrapper, true> in1,
		PositIntermediateFormat<N, WES, Wrapper, true> in2
){
	constexpr auto S_WF = PositDim<N, WES>::WF;
	constexpr auto S_WE = PositDim<N, WES>::WE;
	constexpr auto S_WES = WES;
	constexpr auto K_SIZE = S_WE - S_WES;

	//Sort in order to have exponent of in1 greater than exponent of in2
	auto in2IsZero = in2.isZero();
	auto in1IsZero = in1.isZero();
	auto oneIsZero = in1IsZero | in2IsZero;

	auto input2Significand = in2.getSignedSignificand();
	auto input1Significand = in1.getSignedSignificand();

	auto exp1 = in1.getExp().as_signed();
	auto exp2 = in2.getExp().as_signed();

	auto in1IsGreater = (in1IsZero.invert() & (exp1 > exp2)) | in2IsZero;

	auto subExpOp1 = Wrapper<S_WE, true>::mux(in1IsGreater, exp1, exp2);
	auto subExpOp2 = Wrapper<S_WE, true>::mux(in1IsGreater, exp2, exp1);
	auto mostSignificantSignificand = Wrapper<S_WF+2, false>::mux(in1IsGreater, input1Significand, input2Significand);
	auto lessSignificantSignificand = Wrapper<S_WF+2, false>::mux(in1IsGreater, input2Significand, input1Significand);

	auto mostSignifSign = Wrapper<1, false>::mux(in1IsGreater, in1.getSignBit(), in2.getSignBit());
	auto lessSignifSign = Wrapper<1, false>::mux(in1IsGreater, in2.getSignBit(), in1.getSignBit());

	// Relative shift of exponents

	auto shiftValue = subExpOp1.modularSub(subExpOp2).as_unsigned();
	auto shiftedSignificand = shifter_sticky(
				lessSignificantSignificand.concatenate(Wrapper<2, false>{0}),
				shiftValue,
				lessSignifSign
			);

	Wrapper<S_WF + 2, false> shiftedTop = shiftedSignificand.template slice<S_WF+2+2, 3>();
	auto guards = shiftedSignificand.template slice<2, 1>();
	auto sticky_low = shiftedSignificand.template get<0>();

	auto addOp1 = mostSignifSign.concatenate(mostSignificantSignificand);
	auto addOp2 = lessSignifSign.concatenate(shiftedTop);

	auto addRes = addOp1.modularAdd(addOp2);
	auto toCount = addRes.template get<S_WF+2>();
	auto usefulRes = addRes.template slice<S_WF+1, 0>();

	auto lzoc_shifted = LZOC_shift<S_WF+4, S_WF+4>(usefulRes.concatenate(guards), toCount);


	auto lzoc = lzoc_shifted.lzoc;
	auto frac = lzoc_shifted.shifted.template slice<S_WF+3, 3>();
	auto round = lzoc_shifted.shifted.template get<2>();
	auto sticky = sticky_low.bitwise_or(lzoc_shifted.shifted.template get<1>().bitwise_or(lzoc_shifted.shifted.template get<0>()));

	auto is_zero = toCount.invert().bitwise_and(lzoc == Wrapper<Static_Val<S_WF+4>::_storage, false>{S_WF+4});
	auto non_zero_exp = subExpOp1.subWithCarry(lzoc.template leftpad<S_WE>().as_signed(), {1}).template slice<S_WE-1, 0>();
	auto final_exp = Wrapper<S_WE, false>::mux(
					is_zero,
					{0},
					non_zero_exp
				).as_unsigned();

	auto isResultNar = in1.getIsNaR().bitwise_or(in2.getIsNaR());

	PositIntermediateFormat<N, WES, Wrapper, false> result {
				round,
				sticky,
				isResultNar,
				final_exp,
				toCount,
				frac.template get<S_WF>(),
				frac.template slice<S_WF-1, 0>()
	};
#ifdef POSIT_ADDER_DEBUG
	cerr << "=== POSIT_ADD ===" << endl;
	cerr << "in1IsZero: " << to_string(in1IsZero) << endl;
	cerr << "in2IsZero: " << to_string(in2IsZero) << endl;
	cerr << "oneIsZero: " << to_string(oneIsZero) << endl;
	cerr << "input2Significand: " << to_string(input2Significand) << endl;
	cerr << "input1Significand: " << to_string(input1Significand) << endl;
	cerr << "exp1: " << to_string(exp1) << endl;
	cerr << "exp2: " << to_string(exp2) << endl;
	cerr << "in1IsGreater: " << to_string(in1IsGreater) << endl;
	cerr << "subExpOp1: " << to_string(subExpOp1) << endl;
	cerr << "subExpOp2: " << to_string(subExpOp2) << endl;
	cerr << "mostSignificantSignificand: " << to_string(mostSignificantSignificand) << endl;
	cerr << "lessSignificantSignificand: " << to_string(lessSignificantSignificand) << endl;
	cerr << "mostSignifSign: " << to_string(mostSignifSign) << endl;
	cerr << "lessSignifSign: " << to_string(lessSignifSign) << endl;
	cerr << "shiftValue: " << to_string(shiftValue) << endl;
	cerr << "shiftedSignificand: " << to_string(shiftedSignificand) << endl;
	cerr << "shiftedTop: " << to_string(shiftedTop) << endl;
	cerr << "guards: " << to_string(guards) << endl;
	cerr << "sticky_low: " << to_string(sticky_low) << endl;
	cerr << "addOp1: " << to_string(addOp1) << endl;
	cerr << "addOp2: " << to_string(addOp2) << endl;
	cerr << "addRes: " << to_string(addRes) << endl;
	cerr << "toCount: " << to_string(toCount) << endl;
	cerr << "usefulRes: " << to_string(usefulRes) << endl;
	cerr << "lzoc: " << to_string(lzoc_shifted.lzoc) << endl;
	cerr << "shifted: " << to_string(lzoc_shifted.shifted) << endl;
	cerr << "lzoc: " << to_string(lzoc) << endl;
	cerr << "frac: " << to_string(frac) << endl;
	cerr << "round: " << to_string(round) << endl;
	cerr << "sticky: " << to_string(sticky) << endl;
	cerr << "is_zero: " << to_string(is_zero) << endl;
	cerr << "non_zero_exp: " << to_string(non_zero_exp) << endl;
	cerr << "final_exp: " << to_string(final_exp) << endl;
	cerr << "isResultNar: " << to_string(isResultNar) << endl;
	cerr << "==================" << endl;
#endif
	return result;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, true> posit_add_in_place(
		PositIntermediateFormat<N, WES, Wrapper, true> in1,
		PositIntermediateFormat<N, WES, Wrapper, true> in2
){
	auto interm = posit_add(in1, in2);
	auto res = in_place_rounder(interm);
	return res;
}
