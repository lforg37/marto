#ifndef POSIT_IN_PLACE_ROUND_HPP
#define POSIT_IN_PLACE_ROUND_HPP

#include "posit_dim.hpp"
#include "primitives/shifter_sticky.hpp"
#include "primitives/shifter.hpp"

#ifdef POSIT_ROUNDER_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cout;
#endif

using namespace hint;
using namespace std;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, true>
in_place_rounder(PositIntermediateFormat<N, WES, Wrapper, false> to_round)
{
	constexpr auto S_WF = PositDim<N, WES>::WF;
	constexpr auto S_WE = PositDim<N, WES>::WE;
	constexpr auto S_WES = WES;
	constexpr auto K_SIZE = S_WE - S_WES;

	auto sign = to_round.getSignBit();
	auto current_frac_wo_implicit_bit = to_round.getSignificand().template slice<S_WF-1, 0>();
	auto current_implicit_bit = to_round.getImplicitBit();

	auto exp = to_round.getExp();
	auto k = exp.template slice<S_WE-1, S_WES>(); // Take the exponent bits which should be encoded by the range
	auto low_k = k.template slice<K_SIZE-2, 0>();  //Ignore K sign bit
	auto k_sign = k.template get<K_SIZE - 1>();
	auto absK = Wrapper<K_SIZE-1, false>::mux(k_sign, low_k.invert(), low_k);

	auto smask = Wrapper<WES, false>::generateSequence(sign);
	auto es = exp.template slice<S_WES-1, 0>() ^ smask;

	auto isNeg = k_sign ^ sign;
	auto reverse_and_es = isNeg.concatenate(es);

	auto exp_and_frac = reverse_and_es.concatenate(current_frac_wo_implicit_bit);


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

	auto must_add_1 = guard.bitwise_and((sticky).bitwise_or(round));
	auto rounded_frac = masked_fraction.addWithCarry(mask_round, Wrapper<1, false>{0});

	auto final_frac = Wrapper<S_WF, false>::mux(must_add_1,
								rounded_frac.template slice<S_WF-1,0>(),
								masked_fraction.template slice<S_WF - 1, 0>()
								);

	auto final_es = Wrapper<S_WES, false>::mux(must_add_1,
						rounded_frac.template slice<S_WF+S_WES-1,S_WF>(),
						masked_fraction.template slice<S_WF+S_WES-1,S_WF>()
					) ^ smask;

	auto es_overflow = must_add_1 & (rounded_frac.template get<S_WES + S_WF>() ^ masked_fraction.template get<S_WES + S_WF>());

	auto to_add_high_bit = sign & es_overflow;

	auto to_add_k = Wrapper<K_SIZE - 1, false>::generateSequence(to_add_high_bit).concatenate(es_overflow);
	auto final_k = k.modularAdd(to_add_k);

	auto final_exp = final_k.concatenate(final_es);

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
	cerr << "current_frac_wo_implicit_bit: " << to_string(current_frac_wo_implicit_bit) << endl;
	cerr << "current_implicit_bit: " << to_string(current_implicit_bit) << endl;
	cerr << "exp: " << to_string(exp) << endl;
	cerr << "k: " << to_string(k) << endl;
	cerr << "low_k: " << to_string(low_k) << endl;
	cerr << "k_sign: " << to_string(k_sign) << endl;
	cerr << "absK: " << to_string(absK) << endl;
	cerr << "es: " << to_string(es) << endl;
	cerr << "isNeg: " << to_string(isNeg) << endl;
	cerr << "reverse_and_es: " << to_string(reverse_and_es) << endl;
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
	cerr << "final_es: " << to_string(final_es) << endl;
	cerr << "es_overflow: " << to_string(es_overflow) << endl;
	cerr << "to_add_high_bit: " << to_string(to_add_high_bit) << endl;
	cerr << "to_add_k: " << to_string(to_add_k) << endl;
	cerr << "final_k: " << to_string(final_k) << endl;
	cerr << "final_exp: " << to_string(final_exp) << endl;
	cerr << "is_zero: " << to_string(is_zero) << endl;
	cerr << "final_exp_if_zero: " << to_string(final_exp_if_zero) << endl;
	cerr << "isResultNar: " << to_string(isResultNar) << endl;
	cerr << "=============================" << endl;
#endif

	return result;
}
#endif // POSIT_IN_PLACE_ROUND_HPP
