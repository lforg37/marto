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
inline IEEENumber<WE, WF, Wrapper> ieee_product(IEEENumber<WE, WF, Wrapper> i0, IEEENumber<WE, WF, Wrapper> i1)
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
#endif


	auto i0_exp = i0.getExponent();
	auto i0_exp_adjusted = i0_exp.template slice<WE-1,1>().concatenate(i0_exp.template get<0>().bitwise_or(i0_is_subnormal));
	auto i1_exp = i1.getExponent();
	auto i1_exp_adjusted = i1_exp.template slice<WE-1,1>().concatenate(i1_exp.template get<0>().bitwise_or(i1_is_subnormal));

	auto i0_imp_significand = i0.getLeadBitVal().concatenate(i0.getFractionnalPart());
	auto i1_imp_significand = i1.getLeadBitVal().concatenate(i1.getFractionnalPart());

	auto sign = i0.getSign() xor i1.getSign();
	Wrapper<WE+1, false> exp_sum = i0_exp_adjusted.addWithCarry(i1_exp_adjusted, Wrapper<1, false>{0});
	Wrapper<2*WF+2, false> significand_prod = (i0_imp_significand) * (i1_imp_significand);

	Wrapper<WE+2, false> FP_BIAS{(1<<(WE-1))-1};

#ifdef DEBUG
	cerr << "i0 sign: " << to_string(i0.getSign()) << endl;
	cerr << "i1 sign: " << to_string(i1.getSign()) << endl;
	cerr << "i0 exp: " << to_string(i0_exp) << endl;
	cerr << "i1 exp: " << to_string(i1_exp) << endl;
	cerr << "i0 mantissa: " << to_string(i0_imp_significand) << endl;
	cerr << "i1 mantissa: " << to_string(i1_imp_significand) << endl;
	cerr << "FP_BIAS " << to_string(FP_BIAS) << endl;
	cerr << "sign: " << to_string(sign) << endl;
	cerr << "exp_sum: " << to_string(exp_sum) << endl;
	cerr << "significand_prod " << to_string(significand_prod) << endl;
#endif

	auto entry_to_lzc = Wrapper<WF+1, false>::mux(i0_is_subnormal,
											i0_imp_significand,
											i1_imp_significand
											);
	auto lzc_subnormal = lzoc_wrapper(entry_to_lzc, Wrapper<1, false>{0});

#ifdef DEBUG
	cerr << "to_lzc: " << to_string(entry_to_lzc) << endl;
	cerr << "lzc_subnormal: " << to_string(lzc_subnormal) << endl;
#endif

	Wrapper<WE+2, false> ext_exp_sum = Wrapper<1, false>{0}.concatenate(exp_sum);

	auto unbiased_exp = ext_exp_sum.modularSub(FP_BIAS);
	auto unbiased_is_neg = unbiased_exp.template get<WE+2-1>();
	// (
	// 					Wrapper<2, false>{0}.concatenate(
	// 											Wrapper<WE-1, false>::generateSequence({1})
	// 										).concatenate(
	// 											Wrapper<1, false>::generateSequence({0})
	// 										)
	// 				);



	auto right_shift_value = Wrapper<WE+2, false>{1}.modularSub(unbiased_exp);



	auto unbiased_minus_1 = unbiased_exp.modularSub(Wrapper<WE+2,false>{1});
	auto unbiased_minus_1_red = unbiased_minus_1.template slice<WE-1,0>();
	// auto early_of = unbiased_minus_1.template slice<WE+2-1,WE>().or_reduction();
	auto unbiased_minus_1_is_neg = unbiased_minus_1.template get<WE+2-1>();
	auto ext_lzc = lzc_subnormal.template leftpad<WE>();
	auto left_shift_value = Wrapper<WE, false>::mux(ext_lzc>unbiased_minus_1_red,
												unbiased_minus_1_red,
												ext_lzc);


	auto final_exp = Wrapper<WE, false>::mux(unbiased_minus_1_is_neg,
												Wrapper<WE, false>{0},
												unbiased_minus_1_red
											);
#ifdef DEBUG
	cerr << "unbiased_minus_1 " << to_string(unbiased_minus_1) << endl;
	cerr << "ext_lzc " << to_string(ext_lzc) << endl;
	cerr << "final_exp " << to_string(final_exp) << endl;
	cerr << "left_shift_value " << to_string(left_shift_value) << endl;
	cerr << "unbiased_exp " << to_string(unbiased_exp) << endl;
#endif

	auto must_right_shift = unbiased_minus_1.template get<WE+2-1>();
	auto must_left_shift = must_right_shift.invert();
	auto left_shift = Wrapper<WE+2, false>::mux(must_left_shift,
		Wrapper<WE+2, false>{0}.modularSub(left_shift_value.template leftpad<WE+2>()),
		Wrapper<WE+2, false>{0}
		);

	// auto must_right_shift = must_left_shift.invert();
	// unbiased_is_neg.bitwise_or(unbiased_exp.or_reduction().invert());
	auto right_shift = Wrapper<WE+2, false>::mux(must_right_shift,
		right_shift_value,
		left_shift
		);

#ifdef DEBUG
	cerr << "must_right_shift: \t\t\t" << to_string(must_right_shift) << endl;
	cerr << "must_left_shift: \t\t\t" << to_string(must_left_shift) << endl;
	cerr << "right_shift_value: \t\t\t" << to_string(right_shift_value) << endl;
	cerr << "left_shift: \t\t\t" << to_string(left_shift) << endl;
	cerr << "to sub to : \t\t\t" << to_string(Wrapper<WE+2, false>{2*WF+3-(WF+3)}) << endl;
#endif

	auto shift_value = Wrapper<WE+2, false>{2*WF+3-(WF+3)}.modularAdd(right_shift);
	auto right_shifted_with_sticky = hint::shifter_sticky(significand_prod.concatenate(Wrapper<1, false>{0}), shift_value, Wrapper<1, false>{0});
	auto right_shifted = right_shifted_with_sticky.template slice<2*WF+2+1-1,1>();
	auto sticky_from_shift = right_shifted_with_sticky.template get<0>();

#ifdef DEBUG
	cerr << "shift_input: \t\t\t" << to_string(significand_prod.concatenate(Wrapper<1, false>{0})) << endl;
	cerr << "shift_value: \t\t\t" << to_string(shift_value) << endl;
	cerr << "right_shifted_with_sticky: \t" << to_string(right_shifted_with_sticky) << endl;
	cerr << "right_shifted: \t\t\t" << to_string(right_shifted) << endl;
	cerr << "sticky_from_shift: \t\t\t" << to_string(sticky_from_shift) << endl;
#endif

	auto needs_extra_shift = right_shifted.template get<WF+2>();
	auto lead_exact_prod = Wrapper<1,false>::mux(needs_extra_shift,
												needs_extra_shift,
												right_shifted.template get<WF+2-1>()
											);

	auto final_significand = Wrapper<WF,false>::mux(needs_extra_shift,
												right_shifted.template slice<WF+1,2>(),
												right_shifted.template slice<WF,1>()
											);

	auto round_bit = Wrapper<1,false>::mux(needs_extra_shift,
												right_shifted.template get<1>(),
												right_shifted.template get<0>()
											);

	auto lsb = Wrapper<1,false>::mux(needs_extra_shift,
												right_shifted.template get<2>(),
												right_shifted.template get<1>()
											);

	auto sticky = Wrapper<1,false>::mux(needs_extra_shift,
												sticky_from_shift.bitwise_or(right_shifted.template get<0>()),
												sticky_from_shift
											);

	auto to_sum_for_rounding = round_bit.bitwise_and(sticky).bitwise_or(round_bit.bitwise_and(sticky.invert()).bitwise_and(lsb));
	auto scaled_exp = final_exp.addWithCarry(Wrapper<WE, false>{0}, Wrapper<1, false>{1});
	auto adjusted_exp_if_sub = Wrapper<WE+1, false>::mux(i0_is_subnormal.bitwise_or(i1_is_subnormal).bitwise_and(must_left_shift),
														scaled_exp.modularSub(left_shift_value.template leftpad<WE+1>()),
														scaled_exp);
	auto adjusted_exp_if_rebecomes_normal = Wrapper<WE+1, false>::mux(adjusted_exp_if_sub.template slice<WE,1>().or_reduction().invert().bitwise_and(adjusted_exp_if_sub.template get<0>()).bitwise_and(lead_exact_prod.invert()),
														Wrapper<WE+1, false>{0},
														adjusted_exp_if_sub
														);

	auto adjusted_exp_extra_shift = Wrapper<WE+2, false>::mux(needs_extra_shift,
														adjusted_exp_if_rebecomes_normal.addWithCarry(Wrapper<WE+1, false>{0}, Wrapper<1, false>{1}),
														Wrapper<1, false>{0}.concatenate(adjusted_exp_if_rebecomes_normal)
														);
	auto threshold = Wrapper<WE+2, false>{(1<<WE)-1};
	auto overflow = (adjusted_exp_extra_shift>=threshold).bitwise_or(unbiased_minus_1.template get<WE+1>().invert().bitwise_and(unbiased_minus_1.template get<WE>()));
#ifdef DEBUG
	cerr << "scaled_exp " << to_string(scaled_exp) << endl;
	cerr << "threshold " << to_string(threshold) << endl;
	cerr << "overflow " << to_string(overflow) << endl;
#endif
#ifdef DEBUG

	cerr << "lead_exact_prod: \t\t" << to_string(lead_exact_prod) << endl;
	cerr << "round_bit: \t\t\t" << to_string(round_bit) << endl;
	cerr << "lsb: \t\t\t\t" << to_string(lsb) << endl;
	cerr << "to_sum_for_rounding: \t\t" << to_string(to_sum_for_rounding) << endl;
#endif



	auto rounded_result = sign.concatenate(adjusted_exp_extra_shift.template slice<WE-1,0>()).concatenate(final_significand).modularAdd(to_sum_for_rounding.template leftpad<WE+WF+1>());

#ifdef DEBUG
	cerr << "rounded_result: \t " << to_string(rounded_result) << endl;
	cerr << "adjusted_exp_if_sub: \t " << to_string(adjusted_exp_if_sub) << endl;
	cerr << "adjusted_exp_if_rebecomes_normal: \t " << to_string(adjusted_exp_if_rebecomes_normal) << endl;
	cerr << "adjusted_exp_extra_shift: \t " << to_string(adjusted_exp_extra_shift) << endl;
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

	auto result_if_zero = Wrapper<1+WE+WF, false>::mux(i0_is_zero.bitwise_or(i1_is_zero),
														sign.concatenate(Wrapper<WE+WF, false>{0}),
														rounded_result
														);
	auto result_if_NaN = Wrapper<1+WE+WF, false>::mux(result_is_NaN,
														NaN,
														result_if_zero
														);
	auto result_if_infty = Wrapper<1+WE+WF, false>::mux(result_is_infty,
													signed_infty,
													result_if_NaN
													);
	return result_if_infty;

}
#endif // IEEE_MULTIPLIER_HPP
