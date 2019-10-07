#ifndef IEEE_MULTIPLIER_HPP
#define IEEE_MULTIPLIER_HPP

#include "ieeetype.hpp"

#include "primitives/lzoc.hpp"
#include "primitives/shifter_sticky.hpp"

// #define DEBUG
#ifdef DEBUG
	#include "tools/printing.hpp"
#endif

using namespace hint;

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
IEEENumber<WE, WF, Wrapper> ieee_product(IEEENumber<WE, WF, Wrapper> i0, IEEENumber<WE, WF, Wrapper> i1)
{
	auto i0_frac_is_zero = i0.getFractionnalPart().or_reduction().invert();
	auto i1_frac_is_zero = i1.getFractionnalPart().or_reduction().invert();
	auto i0_exp_is_zero = i0.getExponent().or_reduction().invert();
	auto i1_exp_is_zero = i1.getExponent().or_reduction().invert();
	auto i0_is_zero = i0_exp_is_zero.bitwise_and(i0_frac_is_zero);
	auto i1_is_zero = i1_exp_is_zero.bitwise_and(i1_frac_is_zero);
	auto i0_is_subnormal = i0_exp_is_zero;
	auto i1_is_subnormal = i1_exp_is_zero;
	auto i0_is_infty = i0.getExponent().and_reduction().bitwise_and(i0_frac_is_zero);
	auto i1_is_infty = i1.getExponent().and_reduction().bitwise_and(i1_frac_is_zero);




#ifdef DEBUG
	cerr << "i0: " << to_string(i0) << endl;
	cerr << "i1: " << to_string(i1) << endl;
	cerr << "i0 sign: " << to_string(i0.getSign()) << endl;
	cerr << "i1 sign: " << to_string(i1.getSign()) << endl;
#endif

	auto i0_imp_significand = i0.getLeadBitVal().concatenate(i0.getFractionnalPart());
	auto i1_imp_significand = i1.getLeadBitVal().concatenate(i1.getFractionnalPart());

#ifdef DEBUG
	cerr << "i0 mantissa: " << to_string(i0_imp_significand) << endl;
	cerr << "i1 mantissa: " << to_string(i1_imp_significand) << endl;
#endif

	auto entry_to_lzc = Wrapper<WF+1, false>::mux(i0_is_subnormal,
											i0_imp_significand,
											i1_imp_significand
											);
	auto lzc_subnormal = lzoc_wrapper(entry_to_lzc, Wrapper<1, false>{0});

#ifdef DEBUG
	cerr << "lzc_subnormal: " << to_string(lzc_subnormal) << endl;
#endif



	auto sign = i0.getSign() xor i1.getSign();

#ifdef DEBUG
	cerr << "sign: " << to_string(sign) << endl;
#endif

	Wrapper<WE+2, false> FP_BIAS{(1<<(WE-1))-1};

#ifdef DEBUG
	cerr << "FP_BIAS " << to_string(FP_BIAS) << endl;
	cerr << "i0 exp: \t " << to_string(i0.getExponent()) << endl;
	cerr << "i1 exp: \t " << to_string(i1.getExponent()) << endl;
#endif

	Wrapper<WE+1, false> exp_sum = i0.getExponent().addWithCarry(i1.getExponent(), Wrapper<1, false>{0});
	Wrapper<WE+1, false> exp_sum_plus_one = i0.getExponent().addWithCarry(i1.getExponent(), Wrapper<1, false>{1});
	Wrapper<2*WF+2, false> significand_prod = (i0_imp_significand) * (i1_imp_significand);

	Wrapper<WE+1, false> exp_sum_to_take = Wrapper<WE+1, false>::mux(significand_prod.template get<2*WF+2-1>(), exp_sum_plus_one, exp_sum);
#ifdef DEBUG
	cerr << "exp_sum_to_take " << to_string(exp_sum_to_take) << endl;
#endif
	Wrapper<WE+2, false> ext_exp_sum = Wrapper<1, false>{0}.concatenate(exp_sum_to_take);

	auto unbiased_exp = ext_exp_sum.modularSub(FP_BIAS);
	auto unbiased_is_neg = unbiased_exp.template get<WE+2-1>();
	auto threshold = (Wrapper<2, false>{0}.concatenate(Wrapper<WE, false>::generateSequence({1})));
	auto overflow = (unbiased_exp>=threshold).bitwise_and(unbiased_exp.template get<WE+2-1>().invert());
#ifdef DEBUG
	cerr << "threshold " << to_string(threshold) << endl;
	cerr << "overflow " << to_string(overflow) << endl;
#endif

	auto to_add_to_right_shift = significand_prod.template get<2*WF+2-1>().concatenate(
							significand_prod.template get<2*WF+2-2>().bitwise_and(significand_prod.template get<2*WF+2-1>().invert()).bitwise_and(i0_is_subnormal.bitwise_or(i1_is_subnormal).invert())
							);
#ifdef DEBUG
	cerr << "to_add_to_right_shift " << to_string(to_add_to_right_shift) << endl;
#endif


	auto right_shift_value = Wrapper<WE, false>{0}.concatenate(to_add_to_right_shift).modularSub(unbiased_exp);



	auto final_exp = Wrapper<WE, false>::mux(unbiased_is_neg,
												Wrapper<WE, false>{0},
												unbiased_exp.template slice<WE-1,0>()
											);

	auto ext_lzc = lzc_subnormal.template leftpad<WE>();
	auto left_shift_value = Wrapper<WE, false>::mux(ext_lzc>final_exp,
												final_exp,
												ext_lzc);

#ifdef DEBUG
	cerr << "unbiased_exp " << to_string(unbiased_exp) << endl;
	cerr << "right_shift_value " << to_string(right_shift_value) << endl;
	cerr << "left_shift_value " << to_string(left_shift_value) << endl;
#endif

	auto adjusted_exp_if_sub = Wrapper<WE, false>::mux(i0_is_subnormal.bitwise_or(i1_is_subnormal).bitwise_and(unbiased_is_neg.invert()),
														final_exp.modularSub(left_shift_value),
														final_exp);
	auto adjusted_exp_if_zero = Wrapper<WE, false>::mux(i0_is_zero.bitwise_or(i1_is_zero),
														Wrapper<WE, false>{0},
														adjusted_exp_if_sub);

	auto must_left_shift = i0_is_subnormal.bitwise_or(i1_is_subnormal);
	auto left_shift = Wrapper<WE+2, false>::mux(must_left_shift,
		Wrapper<WE+2, false>{0}.modularSub(left_shift_value.template leftpad<WE+2>()),
		Wrapper<WE+2, false>{0}
		);

	auto must_right_shift = unbiased_is_neg.bitwise_or(unbiased_exp.or_reduction().invert());
	auto right_shift = Wrapper<WE+2, false>::mux(must_right_shift,
		right_shift_value,
		left_shift
		);

	auto shift_value = Wrapper<WE+2, false>{2*WF+3-(WF+3)}.modularAdd(right_shift);
	auto right_shifted_with_sticky = hint::shifter_sticky(significand_prod.concatenate(Wrapper<1, false>{0}), shift_value, Wrapper<1, false>{0});
	auto right_shifted = right_shifted_with_sticky.template slice<2*WF+2+1-1,1>();
	auto sticky_from_shift = right_shifted_with_sticky.template get<0>();

	auto lead_exact_prod = right_shifted.template get<WF+2>();
	auto second_lead_exact_prod = right_shifted.template get<WF+2-1>();
	auto final_significand = Wrapper<WF, false>::mux(lead_exact_prod,
														right_shifted.template slice<WF+1,2>(),
														right_shifted.template slice<WF,1>()
														);
	auto extra_sticky_bit = Wrapper<1, false>::mux(lead_exact_prod,
														right_shifted.template get<0>(),
														Wrapper<1, false>{0}
														);
	auto final_sticky = sticky_from_shift.bitwise_or(extra_sticky_bit);
	auto round_bit = Wrapper<1, false>::mux(lead_exact_prod,
														right_shifted.template get<1>(),
														right_shifted.template get<0>()
														);
	auto lsb = Wrapper<1, false>::mux(lead_exact_prod,
													right_shifted.template get<2>(),
													right_shifted.template get<1>()
													);
	auto to_sum_for_rounding = round_bit.bitwise_and(final_sticky).bitwise_or(round_bit.bitwise_and(final_sticky.invert()).bitwise_and(lsb));


#ifdef DEBUG
	cerr << "final_sticky: \t\t\t" << to_string(final_sticky) << endl;
	cerr << "round_bit: \t\t\t" << to_string(round_bit) << endl;
	cerr << "lsb: \t\t\t\t" << to_string(lsb) << endl;
	cerr << "to_sum_for_rounding: \t\t" << to_string(to_sum_for_rounding) << endl;
#endif



#ifdef DEBUG
	cerr << "base: \t\t\t" << to_string(Wrapper<WE+2, false>{2*WF+2-(WF+3)}) << endl;
	cerr << "add1: \t\t\t" << to_string(significand_prod.template get<2*WF+2-1>()) << endl;
	cerr << "left: \t\t\t" << to_string(left_shift) << endl;
	cerr << "right: \t\t\t" << to_string(right_shift) << endl;
	cerr << "shift_value: \t\t\t" << to_string(shift_value) << endl;
	cerr << "right_shifted: \t\t\t" << to_string(right_shifted) << endl;
	cerr << "final_significand: \t\t\t" << to_string(final_significand) << endl;
#endif


	auto adjusted_exp_if_rebecomes_normal = Wrapper<WE, false>::mux(i0_is_subnormal.bitwise_or(i1_is_subnormal),
														adjusted_exp_if_zero.modularAdd(Wrapper<WE-2, false>::generateSequence({0}).concatenate(lead_exact_prod).concatenate(lead_exact_prod.invert().bitwise_and(second_lead_exact_prod))),
														adjusted_exp_if_zero);

	auto rounded_result = sign.concatenate(adjusted_exp_if_rebecomes_normal).concatenate(final_significand).modularAdd(to_sum_for_rounding.template leftpad<WE+WF+1>());

#ifdef DEBUG
	cerr << "bias: \t\t" << to_string(FP_BIAS) << endl;
	cerr << "unbiased_exp: \t" << to_string(unbiased_exp) << endl;
	cerr << "final_exp: \t " << to_string(final_exp) << endl;
	cerr << "adjusted_exp_if_sub: \t " << to_string(adjusted_exp_if_sub) << endl;
	cerr << "adjusted_exp_if_zero: \t " << to_string(adjusted_exp_if_zero) << endl;
	cerr << "adjusted_exp_if_rebecomes_normal: \t " << to_string(adjusted_exp_if_rebecomes_normal) << endl;
	cerr << "significand prod: " << to_string(significand_prod) << endl;
	cerr << "final_significand: " << to_string(final_significand) << endl;
#endif


	auto result_is_NaN = i0.isNaN().bitwise_or(i1.isNaN()).bitwise_or(i0_is_zero.bitwise_and(i1_is_infty)).bitwise_or(i1_is_zero.bitwise_and(i0_is_infty));
	auto NaN = Wrapper<1, false>{0}.concatenate(Wrapper<WE+WF, false>::generateSequence({1}));
	auto result_is_infty = i0_is_infty.bitwise_or(i1_is_infty).bitwise_and(result_is_NaN.invert()).bitwise_or(overflow.bitwise_and(result_is_NaN.invert()));
	auto signed_infty = Wrapper<1, false>{sign}.concatenate(Wrapper<WE, false>::generateSequence({1})).concatenate(Wrapper<WF, false>::generateSequence({0}));

#ifdef DEBUG
	cerr << "result_is_NaN: \t\t" << to_string(result_is_NaN) << endl;
	cerr << "result_is_infty: \t" << to_string(result_is_infty) << endl;
#endif

	auto result_if_NaN = Wrapper<1+WE+WF, false>::mux(result_is_NaN,
														NaN,
														rounded_result
														);
	auto result_if_infty = Wrapper<1+WE+WF, false>::mux(result_is_infty,
													signed_infty,
													result_if_NaN
													);
	return result_if_infty;

}
#endif // IEEE_MULTIPLIER_HPP
