#pragma once
#include <cstdio>
#include "tools/static_math.hpp"
#include "primitives/shifter_sticky.hpp"
#include "primitives/lzoc_shifter.hpp"
#include "primitives/table.hpp"
#include "helpers/bit_sequence_generator.hpp"

#include "posit_dim.hpp"

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

	auto min_is_zero_mask = Wrapper<S_WE, false>::generateSequence(oneIsZero.invert());
	auto shiftValue = subExpOp1.modularSub(subExpOp2).as_unsigned() & min_is_zero_mask;
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
	cerr << "lzoc_shifted: " << to_string(lzoc_shifted) << endl;
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
	auto res_sign = addRes.template get<S_WF+2>();
	auto usefulRes = addRes.template slice<S_WF+1, 0>();

	auto lzoc_shifted = LZOC_shift<S_WF+4, S_WF+4>(usefulRes.concatenate(guards), res_sign);


	auto lzoc = lzoc_shifted.template slice<S_WF + 3 + Static_Val<S_WF+4>::_storage, S_WF + 4>();
	auto current_frac = lzoc_shifted.template slice<S_WF+3, 3>();
	auto current_exp = subExpOp1.subWithCarry(lzoc.template leftpad<S_WE>(), {1}).template slice<S_WE-1, 0>();

	// cerr << "old guard: " << to_string(lzoc_shifted.template get<1>()) << endl;
	// cerr << "old sticky: " << to_string(lzoc_shifted.template get<0>()) << endl;


	auto expWoBias = current_exp.modularSub(Wrapper<S_WE, false>{PositDim<N, WES>::EXP_BIAS});
	auto k = expWoBias.template slice<S_WE-1, S_WES>();
	// auto es = expWoBias.template slice<S_WES-1, 0>();
	auto low_k = k.template slice<K_SIZE-2, 0>();
	auto absK = Wrapper<K_SIZE-1, false>::mux(
				k.template get<K_SIZE-1>(),
				low_k.invert(),
				low_k
			);


	Wrapper<1, false> zero_one{1};
	Wrapper<1, false> one_zero{0};


	auto sign_sequence_wes = Wrapper<S_WES, false>::generateSequence(res_sign);
	auto es = expWoBias.template slice<S_WES-1,0>().bitwise_xor(sign_sequence_wes);

	auto isNegative = k.template get<K_SIZE-1>().bitwise_xor(res_sign);

	auto reverse_and_es = Wrapper<S_WES+1, false>::mux(isNegative,
							zero_one.concatenate(es),
							one_zero.concatenate(es)
							);

	auto current_implicit_bit = current_frac.template get<S_WF+1-1>();
	auto current_frac_wo_implicit_bit = current_frac.template slice<S_WF+1-2, 0>();

	auto exp_and_frac = reverse_and_es.concatenate(current_frac_wo_implicit_bit);


	auto mask_top_ones = Wrapper<1, false>::generateSequence({1});
	auto mask_top_zeros = Wrapper<1, false>::generateSequence({0});
	// auto mask_top_zeros_short = Wrapper<2-1, false>::generateSequence({0});

	auto mask_remove_rounded_bits_bottom = hint::one_then_zeros<K_SIZE-1, (1<<(K_SIZE-1)), Wrapper>(absK).template slice<S_WES+S_WF-1, 0>();
	auto mask_remove_rounded_bits = mask_top_ones.concatenate(mask_remove_rounded_bits_bottom);

	auto mask_round_bottom = hint::one_one<K_SIZE-1, (1<<(K_SIZE-1)), (S_WES+S_WF), Wrapper>(absK).template slice<S_WES+S_WF+1-1,0>();
	auto mask_round = mask_round_bottom;
	auto mask_guard = Wrapper<1, false>{0}.concatenate(mask_round.template slice<1+S_WES+S_WF-1,1>());
	auto mask_stickies = (Wrapper<1, false>{1}.concatenate(mask_remove_rounded_bits.template slice<1+S_WES+S_WF-1,1>())).invert();


	auto masked_fraction = exp_and_frac.bitwise_and(mask_remove_rounded_bits);

	auto round_bits = exp_and_frac.bitwise_and(mask_round);
	auto round = round_bits.or_reduction();

	auto guard_bits = exp_and_frac.bitwise_and(mask_guard);
	auto guard_from_mask = guard_bits.or_reduction();


	auto sticky_bits = exp_and_frac.bitwise_and(mask_stickies);
	auto sticky_from_mask = sticky_bits.or_reduction();


	auto guard_is_in_fraction = mask_guard.or_reduction();
	auto guard_outside_fraction = lzoc_shifted.template get<2>();
	auto guard = Wrapper<1, false>::mux(
						guard_is_in_fraction,
						guard_from_mask,
						guard_outside_fraction
					);

	auto sticky_outside_fraction = sticky_low.bitwise_or(lzoc_shifted.template get<1>().bitwise_or(lzoc_shifted.template get<0>()));
	auto sticky = Wrapper<1, false>::mux(
						guard_is_in_fraction,
						sticky_from_mask.bitwise_or(guard_outside_fraction).bitwise_or(sticky_outside_fraction),
						sticky_outside_fraction
					);

	auto must_add_1 = guard.bitwise_and((sticky).bitwise_or(round));

	auto rounded_frac = masked_fraction.addWithCarry(mask_round, Wrapper<1, false>{0});

	auto final_frac = Wrapper<S_WF, false>::mux(must_add_1,
								rounded_frac.template slice<S_WF-1,0>(),
								current_frac_wo_implicit_bit
								);

	auto exp_increased = rounded_frac.template get<S_WF>().bitwise_xor(es.template get<0>());

	auto to_add_to_exp = Wrapper<S_WE, false>::mux(res_sign,
								Wrapper<S_WE, false>::generateSequence(Wrapper<1, false>{1}),
								Wrapper<S_WE, false>{1}
								);
	auto final_exp = Wrapper<S_WE, false>::mux(exp_increased.bitwise_and(must_add_1),
								current_exp.modularAdd(to_add_to_exp),
								//rounded_exp.bitwise_xor(sign_sequence_we),
								current_exp
								);


	auto is_zero = res_sign.invert().bitwise_and(lzoc == Wrapper<Static_Val<S_WF+4>::_storage, false>{S_WF+4});
	auto final_exp_if_zero = Wrapper<S_WE, false>::mux(
					is_zero,
					{0},
					final_exp
				);

	auto isResultNar = in1.getIsNaR().bitwise_or(in2.getIsNaR());

	PositIntermediateFormat<N, WES, Wrapper, true> result {
				isResultNar,
				final_exp_if_zero,
				res_sign,
				current_implicit_bit,
				final_frac.template slice<S_WF-1, 0>()
	};
	return result;
}
