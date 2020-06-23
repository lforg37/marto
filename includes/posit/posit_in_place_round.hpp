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

	auto current_exp = to_round.getExp();
	auto expWoBias = current_exp.modularSub({{PositDim<N, WES>::EXP_BIAS}});
	auto k = expWoBias.template slice<S_WE-1, S_WES>(); // Take the exponent bits which should be encoded by the range
	auto low_k = k.template slice<K_SIZE-2, 0>();  //Ignore K sign bit
	auto k_sign = k.template get<K_SIZE - 1>();
	auto absK = Wrapper<K_SIZE-1, false>::mux(k_sign, low_k.invert(), low_k);

	auto expWoBias_low = expWoBias.template slice<S_WES-1, 0>();
	auto es = Wrapper<S_WES, false>::mux(sign, expWoBias_low.invert(), expWoBias_low);

	auto isNeg = k_sign ^ sign;
	auto reverse_and_es = isNeg.concatenate(es);

	auto exp_and_frac = reverse_and_es.concatenate(current_frac_wo_implicit_bit);


	Wrapper<1, false> mask_top_ones{1}, mask_top_zeros{1};

	auto mask_remove_rounded_bits_bottom = hint::one_then_zeros<K_SIZE-1, (1<<(K_SIZE-1)), Wrapper>(absK).template slice<S_WES+S_WF-1, 0>();
	auto mask_remove_rounded_bits = mask_top_ones.concatenate(mask_remove_rounded_bits_bottom);

	auto mask_round = hint::one_one<K_SIZE-1, (1<<(K_SIZE-1)), (S_WES+S_WF), Wrapper>(absK).template slice<S_WES+S_WF+1-1,0>();
	auto mask_guard = Wrapper<1, false>{0}.concatenate(mask_round.template slice<1+S_WES+S_WF-1,1>());
	auto mask_stickies = (Wrapper<1, false>{1}.concatenate(mask_remove_rounded_bits.template slice<1+S_WES+S_WF-1,1>())).invert();

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

	auto exp_increased = rounded_frac.template get<S_WF>().bitwise_xor(es.template get<0>());

	auto to_add_to_exp = Wrapper<S_WE-1, false>::generateSequence(sign).concatenate(Wrapper<1, false>{1});

	auto final_exp = Wrapper<S_WE, false>::mux(exp_increased.bitwise_and(must_add_1),
								current_exp.modularAdd(to_add_to_exp),
								//rounded_exp.bitwise_xor(sign_sequence_we),
								current_exp
								);

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
	cout << "sign: " << to_string(sign) << endl;
	cout << "current_implicit_bit: " << to_string(current_implicit_bit) << endl;
	cout << "frac_wo_imp_bit: " << to_string(current_frac_wo_implicit_bit) << endl;
	cout << "current_exp: " << to_string(current_exp) << endl;
	cout << "exp_wo_bias: " << to_string(expWoBias) << endl;
	cout << "k: " << to_string(k) << endl;
	cout << "abs_k: " << to_string(absK) << endl;
	cout << "es: " << to_string(es) << endl;
	cout << "isNeg: " << to_string(isNeg) << endl;
	cout << "reverse_and_es: " << to_string(reverse_and_es) << endl;
	cout << "exp_and_frac:" << to_string(exp_and_frac) << endl;
	cout << "mask_remove_rb: " << to_string(mask_remove_rounded_bits) << endl;
	cout << "mask_round: " << to_string(mask_round) << endl;
	cout << "mask_guard: " << to_string(mask_guard) << endl;
	cout << "mask_stickies: " << to_string(mask_stickies) << endl;
	cout << "masked_fraction :" << to_string(masked_fraction) << endl;
	cout << "round: " << to_string(round) << endl;
	cout << "guard_from_mask: " << to_string(guard_from_mask) << endl;
	cout << "sticky_from_mask: " << to_string(sticky_from_mask) << endl;
	cout << "guard_is_in_fraction: " << to_string(guard_is_in_fraction) << endl;
	cout << "guard: " << to_string(guard) << endl;
	cout << "sticky: " << to_string(sticky) << endl;
	cout << "must_add1:" << to_string(must_add_1) << endl;
	cout << "rounded_frac: " << to_string(rounded_frac) << endl;
	cout << "final_frac: " << to_string(final_frac) << endl;

#endif

	return result;
}
#endif // POSIT_IN_PLACE_ROUND_HPP
