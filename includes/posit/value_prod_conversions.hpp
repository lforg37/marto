#pragma once

#include <cstdio>
#include <tools/printing.hpp>

//#include "posit_dim.hpp"

//using hint::to_string;

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


template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, true> PositProd_to_PositIF_in_place_rounding(PositProd<N, WES, Wrapper> val)
{

	constexpr unsigned int S_WF = PositDim<N, WES>::WF;
	constexpr unsigned int S_WE = PositDim<N, WES>::WE;
	constexpr unsigned int S_WES = PositDim<N, WES>::WES;
	constexpr unsigned int PROD_WF = PositDim<N, WES>::ProdSignificandSize;
	constexpr unsigned int PROD_WE = PositDim<N, WES>::ProdExpSize;
	constexpr unsigned int K_SIZE = S_WE - S_WES;

	auto isNaR = val.getIsNaR();
	auto fraction = val.getSignificand();
	auto implicitBit = fraction.template get<PROD_WF-1>();
	auto sign = val.getSignBit();

	auto current_frac = fraction.template slice<PROD_WF - 2, PROD_WF-1-S_WF>();
	//cerr << to_string(current_frac) << endl;

	auto resultGuardBit = fraction.template get<PROD_WF-2-S_WF>();

	auto resultStickyBit = fraction.template slice<PROD_WF-3 -S_WF, 0>().or_reduction();

	auto expt = val.getExp();

	auto expt_no_bias = expt.template leftpad<PROD_WE + 1>().modularSub({PositDim<N, WES>::EXP_BIAS-1});

///// NEW
	auto current_exp = expt_no_bias.template slice<PositDim<N, WES>::WE - 1, 0>();
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


	auto sign_sequence_wes = Wrapper<S_WES, false>::generateSequence(sign);
	auto es = expWoBias.template slice<S_WES-1,0>().bitwise_xor(sign_sequence_wes);

	auto isNegative = k.template get<K_SIZE-1>().bitwise_xor(sign);

	auto reverse_and_es = Wrapper<S_WES+1, false>::mux(isNegative,
							zero_one.concatenate(es),
							one_zero.concatenate(es)
							);

	// auto current_implicit_bit = current_frac.template get<S_WF+1-1>();
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

	// cerr << "exp_and_frac:             " << to_string(exp_and_frac) << endl;
	// cerr << "mask_remove_rounded_bits: " << to_string(mask_remove_rounded_bits) << endl;
	// cerr << "mask_round:               " << to_string(mask_round) << endl;
	// cerr << "mask_guard:               " << to_string(mask_guard) << endl;
	// cerr << "mask_stickies:            " << to_string(mask_stickies) << endl;

	auto masked_fraction = exp_and_frac.bitwise_and(mask_remove_rounded_bits);

	auto round_bits = exp_and_frac.bitwise_and(mask_round);
	auto round = round_bits.or_reduction();

	auto guard_bits = exp_and_frac.bitwise_and(mask_guard);
	auto guard_from_mask = guard_bits.or_reduction();


	auto sticky_bits = exp_and_frac.bitwise_and(mask_stickies);
	auto sticky_from_mask = sticky_bits.or_reduction();


	auto guard_is_in_fraction = mask_guard.or_reduction();
	auto guard_outside_fraction = resultGuardBit;
	auto guard = Wrapper<1, false>::mux(
						guard_is_in_fraction,
						guard_from_mask,
						guard_outside_fraction
					);

	auto sticky_outside_fraction = resultStickyBit;
	auto sticky = Wrapper<1, false>::mux(
						guard_is_in_fraction,
						sticky_from_mask.bitwise_or(guard_outside_fraction).bitwise_or(sticky_outside_fraction),
						sticky_outside_fraction
					);

	auto must_add_1 = guard.bitwise_and((sticky).bitwise_or(round));
	// cerr << "must_add_1: " << to_string(must_add_1) << endl;
	auto rounded_frac = masked_fraction.addWithCarry(mask_round, Wrapper<1, false>{0});

	auto final_frac = Wrapper<S_WF, false>::mux(must_add_1,
								rounded_frac.template slice<S_WF-1,0>(),
								current_frac_wo_implicit_bit
								);

	auto exp_increased = rounded_frac.template get<S_WF>().bitwise_xor(es.template get<0>());

	auto to_add_to_exp = Wrapper<S_WE, false>::mux(sign,
								Wrapper<S_WE, false>::generateSequence(Wrapper<1, false>{1}),
								Wrapper<S_WE, false>{1}
								);
	auto final_exp = Wrapper<S_WE, false>::mux(exp_increased.bitwise_and(must_add_1),
								current_exp.modularAdd(to_add_to_exp),
								//rounded_exp.bitwise_xor(sign_sequence_we),
								current_exp
								);

	// cerr << "final_exp: " << to_string(final_exp) << endl;



//// END NEW







	// hint<1> isZero = not(((hint<4>) val.getSignificand().slice(PositDim<N>::ProdSignificandSize - 1, PositDim<N>::ProdSignificandSize - 4)).or_reduce());
	auto isZero = sign.invert() & implicitBit.invert() & expt.or_reduction().invert();

	auto exp = Wrapper<PROD_WE + 1, false>::mux(
					isZero,
					{0},
					expt_no_bias
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
				PositIntermediateFormat<N, WES, Wrapper, true>(
								isNaR,
								final_exp, //Warning : biased exp
								val.getSignBit(),
								implicitBit,
								final_frac
							)
			);
	/*
	cerr << to_string(resultGuardBit) << endl;
	cerr << to_string(resultStickyBit) << endl;
	cerr << to_string(isNaR) << endl;
	cerr << to_string(exp) << endl;
	cerr << to_string(val.getSignBit()) << endl;
	cerr << to_string(implicitBit) << endl;
	cerr << to_string(current_frac) << endl;
	cerr << to_string(ret) << endl;
	*/
	return ret;



	return PositIntermediateFormat<N, WES, Wrapper, true>{{0}};
}
