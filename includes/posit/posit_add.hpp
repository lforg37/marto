#pragma once
#include <cstdio>
#include "tools/static_math.hpp"
#include "primitives/shifter_sticky.hpp"
#include "primitives/lzoc_shifter.hpp"

#include "posit_dim.hpp"

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
	static constexpr int EXT_SUM_SIZE = Static_Val<S_WF+2 + S_WF +1>::_2pow;
	static constexpr int LOG2_EXT_SUM_SIZE = Static_Val<EXT_SUM_SIZE>::_log2;

	//Sort in order to have exponent of in1 greater than exponent of in2
	auto in1IsGreater = in1.getExp() > in2.getExp();

	auto input2Significand = in2.getSignedSignificand();
	auto input1Significand = in1.getSignedSignificand();

	auto subExpOp1 = Wrapper<S_WE, false>::mux(in1IsGreater, in1.getExp(), in2.getExp());
	auto subExpOp2 = Wrapper<S_WE, false>::mux(in1IsGreater, in2.getExp(), in1.getExp());
	auto mostSignificantSignificand = Wrapper<S_WF+2, false>::mux(in1IsGreater, input1Significand, input2Significand);
	auto lessSignificantSignificand = Wrapper<S_WF+2, false>::mux(in1IsGreater, input2Significand, input1Significand);

	auto mostSignifSign = Wrapper<1, false>::mux(in1IsGreater, in1.getSignBit(), in2.getSignBit());
	auto lessSignifSign = Wrapper<1, false>::mux(in1IsGreater, in2.getSignBit(), in1.getSignBit());

	// Relative shift of exponents
	auto shiftValue = subExpOp1.modularSub(subExpOp2);
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


	auto lzoc = lzoc_shifted.template slice<S_WF + 3 + Static_Val<S_WF+4>::_storage, S_WF + 4>();
	auto frac = lzoc_shifted.template slice<S_WF+3, 3>();
	auto round = lzoc_shifted.template get<2>();
	auto sticky = sticky_low.bitwise_or(lzoc_shifted.template get<1>().bitwise_or(lzoc_shifted.template get<0>()));

	auto is_zero = toCount.invert().bitwise_and(lzoc == Wrapper<Static_Val<S_WF+4>::_storage, false>{S_WF+4});

	auto final_exp = Wrapper<S_WE, false>::mux(
					is_zero,
					{0},
					subExpOp1.subWithCarry(lzoc.template leftpad<S_WE>(), {1}).template slice<S_WE-1, 0>()
				);

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
	return result;
}


template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, true> posit_add_in_place(
		PositIntermediateFormat<N, WES, Wrapper, true> in1,
		PositIntermediateFormat<N, WES, Wrapper, true> in2
){

	// TODO
	cerr << "Use of empty function PositIntermediateFormat<N, WES, Wrapper, true> posit_add_in_place(PositIntermediateFormat<N, WES, Wrapper, true> in1, PositIntermediateFormat<N, WES, Wrapper, true> in2)" << endl;

	return PositIntermediateFormat<N, WES, Wrapper, true>{{0}};
}
