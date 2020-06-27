#ifndef POSIT_IN_PLACE_ROUND_HPP
#define POSIT_IN_PLACE_ROUND_HPP

#include <type_traits>

#include "posit_dim.hpp"
#include "primitives/shifter_sticky.hpp"
#include "primitives/shifter.hpp"
#include "posit_encoder.hpp"

#ifdef POSIT_ROUNDER_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cout;
#endif

using namespace hint;
using namespace std;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<1+WES+PositDim<N, WES>::WF, false> _get_exp_and_frac(
		Wrapper<1, false> valSign,
		Wrapper<1, false> isRangeNeg,
		Wrapper<PositDim<N, WES>::WE, false> exp,
		Wrapper<PositDim<N, WES>::WF, false> frac,
		typename enable_if<PositDim<N, WES>::HAS_ES>::type* = 0
		)
{
	auto smask = Wrapper<WES, false>::generateSequence(valSign);
	auto es = exp.template slice<WES-1, 0>() ^ smask;
	auto reverse_and_es = isRangeNeg.concatenate(es);
	auto exp_and_frac = reverse_and_es.concatenate(frac);
#ifdef POSIT_ROUNDER_DEBUG
	cerr << "=== _get_reverse_and_es (WES != 0) ===" << endl;
	cerr << "smask: " << to_string(smask) << endl;
	cerr << "es: " << to_string(es) << endl;
	cerr << "reverse_and_es: " << to_string(reverse_and_es) << endl;
	cerr << "exp_and_frac: " << to_string(exp_and_frac) << endl;
	cerr << "======================================" << endl;
#endif
	return exp_and_frac;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<1+WES+PositDim<N, WES>::WF, false> _get_exp_and_frac(
		Wrapper<1, false> ,
		Wrapper<1, false> isRangeNeg,
		Wrapper<PositDim<N, WES>::WE, false> ,
		Wrapper<PositDim<N, WES>::WF, false> frac,
		typename enable_if<not PositDim<N, WES>::HAS_ES>::type* = 0
		)
{
	auto exp_and_frac = isRangeNeg.concatenate(frac);
#ifdef POSIT_ROUNDER_DEBUG
	cerr << "=== _get_reverse_and_es (WES == 0) ===" << endl;
	cerr << "exp_and_frac: " << to_string(exp_and_frac) << endl;
	cerr << "======================================" << endl;
#endif
	return exp_and_frac;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<PositDim<N, WES>::WE, false> _get_exp(
			Wrapper<1, false> sign,
			Wrapper<1, false> should_round_up,
			Wrapper<PositDim<N, WES>::WE - WES, false> k,
			Wrapper<1+WES+PositDim<N, WES>::WF, false> truncated_frac,
			Wrapper<1+WES+PositDim<N, WES>::WF, false> ceiled_frac,
			typename enable_if<PositDim<N, WES>::HAS_ES>::type* = 0
		)
{
	constexpr auto WF = PositDim<N, WES>::WF;
	constexpr auto K_SIZE = PositDim<N, WES>::WE - WES;

	auto smask = Wrapper<WES, false>::generateSequence(sign);
	auto final_es = Wrapper<WES, false>::mux(should_round_up,
						ceiled_frac.template slice<WF+WES-1,WF>(),
						truncated_frac.template slice<WF+WES-1,WF>()
					) ^ smask;
	auto es_overflow = should_round_up & (ceiled_frac.template get<WES + WF>() ^ truncated_frac.template get<WES + WF>());

	auto to_add_high_bit = sign & es_overflow;

	auto to_add_k = Wrapper<K_SIZE - 1, false>::generateSequence(to_add_high_bit).concatenate(es_overflow);
	auto final_k = k.modularAdd(to_add_k);

	auto final_exp = final_k.concatenate(final_es);
#ifdef POSIT_ROUNDER_DEBUG
	cerr << "=== _get_final_exp (WES != 0) ===" << endl;
	cerr << "smask: " << to_string(smask) << endl;
	cerr << "final_es: " << to_string(final_es) << endl;
	cerr << "es_overflow: " << to_string(es_overflow) << endl;
	cerr << "to_add_high_bit: " << to_string(to_add_high_bit) << endl;
	cerr << "to_add_k: " << to_string(to_add_k) << endl;
	cerr << "final_k: " << to_string(final_k) << endl;
	cerr << "final_exp: " << to_string(final_exp) << endl;
	cerr << "======================================" << endl;
#endif
	return final_exp;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<PositDim<N, 0>::WE, false> _get_exp(
			Wrapper<1, false> sign,
			Wrapper<1, false> should_round_up,
			Wrapper<PositDim<N, 0>::WE, false> k,
			Wrapper<1+PositDim<N, 0>::WF, false> truncated_frac,
			Wrapper<1+PositDim<N, 0>::WF, false> ceiled_frac,
			typename enable_if<! PositDim<N, WES>::HAS_ES>::type* = 0
		)
{
	constexpr auto WF = PositDim<N, 0>::WF;
	constexpr auto K_SIZE = PositDim<N, 0>::WE;

	auto fracOverflow = should_round_up & (ceiled_frac.template get<WF>() ^ truncated_frac.template get<WF>());

	auto to_add_high_bit = sign & fracOverflow;

	auto to_add_k = Wrapper<K_SIZE - 1, false>::generateSequence(to_add_high_bit).concatenate(fracOverflow);
	auto final_k = k.modularAdd(to_add_k);
#ifdef POSIT_ROUNDER_DEBUG
	cerr << "=== _get_final_exp (WES == 0) ===" << endl;
	cerr << "fracOverflow: " << to_string(fracOverflow) << endl;
	cerr << "to_add_high_bit: " << to_string(to_add_high_bit) << endl;
	cerr << "to_add_k: " << to_string(to_add_k) << endl;
	cerr << "final_k: " << to_string(final_k) << endl;

	cerr << "======================================" << endl;
#endif
	return final_k;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, true>
in_place_rounder(PositIntermediateFormat<N, WES, Wrapper, false> to_round)
{
	constexpr auto S_WF = PositDim<N, WES>::WF;
	constexpr auto S_WE = PositDim<N, WES>::WE;
	constexpr auto S_WES = WES;
	constexpr auto K_SIZE = S_WE - S_WES;

	auto sign = to_round.getSignBit();
	auto current_frac = to_round.getSignificand();
	auto current_frac_wo_implicit_bit = current_frac.template slice<S_WF-1, 0>();
	auto current_implicit_bit = current_frac.template get<S_WF>();

	auto exp = to_round.getExp();
	auto k = exp.template slice<S_WE-1, S_WES>(); // Take the exponent bits which should be encoded by the range
	auto low_k = k.template slice<K_SIZE-2, 0>();  //Ignore K sign bit
	auto k_sign = k.template get<K_SIZE - 1>();
	auto absK = Wrapper<K_SIZE-1, false>::mux(k_sign, low_k.invert(), low_k);

	auto isNeg = k_sign ^ sign;

	auto exp_and_frac = _get_exp_and_frac<N, WES, Wrapper>(sign, isNeg, exp, current_frac_wo_implicit_bit);

	Wrapper<1, false> mask_top_ones{1};

	auto mask_remove_rounded_bits_bottom = hint::one_then_zeros<K_SIZE-1, (1<<(K_SIZE-1)), Wrapper>(absK).template slice<S_WES+S_WF-1, 0>();
	auto mask_remove_rounded_bits = mask_top_ones.concatenate(mask_remove_rounded_bits_bottom);

	auto mask_round = hint::one_one<K_SIZE-1, (1<<(K_SIZE-1)), (S_WES+S_WF), Wrapper>(absK).template slice<S_WES+S_WF,0>();
	auto mask_guard = Wrapper<1, false>{0}.concatenate(mask_round.template slice<S_WES+S_WF,1>());
	auto mask_stickies = (Wrapper<1, false>{1}.concatenate(mask_remove_rounded_bits.template slice<S_WES+S_WF,1>())).invert();

	auto masked_fraction = exp_and_frac.bitwise_and(mask_remove_rounded_bits);

	auto round_bits = exp_and_frac.bitwise_and(mask_round);
	auto round = round_bits.or_reduction();

	auto guard_bits = exp_and_frac.bitwise_and(mask_guard);
	auto guard_from_mask = guard_bits.or_reduction();

	auto sticky_bits = exp_and_frac.bitwise_and(mask_stickies);
	auto sticky_from_mask = sticky_bits.or_reduction();

	//TODO check
	auto guard_is_in_fraction = absK.or_reduction();
	auto guard_outside_fraction = to_round.getGuardBit();
	auto guard = Wrapper<1, false>::mux(
						guard_is_in_fraction,
						guard_from_mask,
						guard_outside_fraction
					);
	auto sticky_outside_fraction = to_round.getStickyBit();
	auto sticky = Wrapper<1, false>::mux(guard_is_in_fraction,
										 guard_outside_fraction | sticky_outside_fraction | sticky_from_mask,
										 sticky_outside_fraction
										 );

	auto roundOverflow = exp_overflow<N, WES, Wrapper>(exp);
	auto forbidRound = isNeg.invert() & roundOverflow;
	auto forceRound = isNeg & roundOverflow;
	auto must_add_1 = (forceRound | guard & (sticky | round)) & forbidRound.invert();
	auto rounded_frac = masked_fraction.modularAdd(mask_round);

	auto final_frac = Wrapper<S_WF, false>::mux(must_add_1,
								rounded_frac.template slice<S_WF-1,0>(),
								masked_fraction.template slice<S_WF - 1, 0>()
								);



	auto final_exp = _get_exp<N, WES, Wrapper>(sign, must_add_1, k, masked_fraction, rounded_frac);

	auto is_zero = to_round.isZero();
	auto final_exp_if_zero = Wrapper<S_WE, false>::mux(
					is_zero,
					{0},
					final_exp
				);

	auto isResultNar = to_round.getIsNaR();

	PositIntermediateFormat<N, WES, Wrapper, true> result {
				isResultNar,
				final_exp_if_zero,
				sign,
				current_implicit_bit,
				final_frac.template slice<S_WF-1, 0>()
	};
#ifdef POSIT_ROUNDER_DEBUG
	cerr << "=== Posit round in place ===" << endl;
	cerr << "Input: " << to_string(to_round) << endl;
	cerr << "sign: " << to_string(sign) << endl;
	cerr << "current_frac: " << to_string(current_frac) << endl;
	cerr << "current_frac_wo_implicit_bit: " << to_string(current_frac_wo_implicit_bit) << endl;
	cerr << "current_implicit_bit: " << to_string(current_implicit_bit) << endl;
	cerr << "exp: " << to_string(exp) << endl;
	cerr << "k: " << to_string(k) << endl;
	cerr << "low_k: " << to_string(low_k) << endl;
	cerr << "k_sign: " << to_string(k_sign) << endl;
	cerr << "absK: " << to_string(absK) << endl;
	cerr << "isNeg: " << to_string(isNeg) << endl;
	cerr << "exp_and_frac: " << to_string(exp_and_frac) << endl;
	cerr << "mask_remove_rounded_bits_bottom: " << to_string(mask_remove_rounded_bits_bottom) << endl;
	cerr << "mask_remove_rounded_bits: " << to_string(mask_remove_rounded_bits) << endl;
	cerr << "mask_round: " << to_string(mask_round) << endl;
	cerr << "mask_guard: " << to_string(mask_guard) << endl;
	cerr << "mask_stickies: " << to_string(mask_stickies) << endl;
	cerr << "masked_fraction: " << to_string(masked_fraction) << endl;
	cerr << "round_bits: " << to_string(round_bits) << endl;
	cerr << "round: " << to_string(round) << endl;
	cerr << "guard_bits: " << to_string(guard_bits) << endl;
	cerr << "guard_from_mask: " << to_string(guard_from_mask) << endl;
	cerr << "sticky_bits: " << to_string(sticky_bits) << endl;
	cerr << "sticky_from_mask: " << to_string(sticky_from_mask) << endl;
	cerr << "guard_is_in_fraction: " << to_string(guard_is_in_fraction) << endl;
	cerr << "guard_outside_fraction: " << to_string(guard_outside_fraction) << endl;
	cerr << "guard: " << to_string(guard) << endl;
	cerr << "sticky_outside_fraction: " << to_string(sticky_outside_fraction) << endl;
	cerr << "sticky: " << to_string(sticky) << endl;
	cerr << "must_add_1: " << to_string(must_add_1) << endl;
	cerr << "rounded_frac: " << to_string(rounded_frac) << endl;
	cerr << "final_frac: " << to_string(final_frac) << endl;
	cerr << "final_exp: " << to_string(final_exp) << endl;
	cerr << "is_zero: " << to_string(is_zero) << endl;
	cerr << "final_exp_if_zero: " << to_string(final_exp_if_zero) << endl;
	cerr << "isResultNar: " << to_string(isResultNar) << endl;
	cerr << "=============================" << endl;
#endif

	return result;
}
#endif // POSIT_IN_PLACE_ROUND_HPP
