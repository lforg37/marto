#pragma once

#include <iostream>

#include "posit_dim.hpp"
#include "primitives/lzoc_shifter.hpp"
//#include "tools/printing.hpp"

//using hint::to_string;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY>
inline PositIntermediateFormat<N, WES, Wrapper, false> quire_to_posit(Quire<N, WES, Wrapper, NB_CARRY> quire)
{
	constexpr int logSize = Static_Val<quire.PositExpRange>::_log2;
	constexpr int allsize = Static_Val<quire.PositExpRange>::_2pow;
	constexpr int padd_width = allsize - quire.PositExpRange;

	auto sign = quire.getSignBit();

	auto upper_low_bits = quire.template slice<
			quire.PositRangeOffset - 1,
			quire.PositRangeOffset - padd_width
		>();

	auto lower_low_bits = quire.template slice<
			quire.PositRangeOffset - padd_width - 1,
			0
		>();

	auto lower_sticky = lower_low_bits.or_reduction();
	auto upper_low_null = upper_low_bits.or_reduction().invert();

	//Are the bits below posit range all null ?
	auto low_bit_is_null = lower_sticky.invert() & upper_low_null;

	auto middle_bits = quire.template slice<
			quire.PositRangeOffset + quire.PositExpRange - 1,
			quire.PositRangeOffset
		>();

	auto middle_s_ext = Wrapper<quire.PositExpRange, false>::generateSequence(sign.invert());
	auto underflow_base = middle_s_ext xor middle_bits;
	auto middle_void_flag = underflow_base.and_reduction();
	auto middle_is_null = middle_bits.or_reduction().invert();

	constexpr unsigned int remainingsize =
		quire.Size - 2 - quire.PositOverflowOffset;

	auto uppercarry = quire.template slice<
			quire.Size - 3,
			quire.PositOverflowOffset
		>();

	auto upper_ext_sign = Wrapper<remainingsize, false>::generateSequence(sign);
	auto overflow = upper_ext_sign xor uppercarry;
	auto high_overflow = overflow.or_reduction();

	auto fin_overflow =
		high_overflow | (sign & middle_is_null & low_bit_is_null);

	auto isZero = sign.invert() &
						fin_overflow.invert() &
						middle_is_null &
						low_bit_is_null;

	auto underflow = high_overflow.invert() & middle_void_flag;

	Wrapper<allsize, false> padded_mid_bits = middle_bits.concatenate(
			upper_low_bits
		);

	auto lzocshifted = hint::LZOC_shift<allsize, allsize-1, false, Wrapper>(padded_mid_bits, sign);
	auto exp = lzocshifted.template slice<logSize + allsize - 1, allsize>();

	auto biased_exp = Wrapper<logSize, false>{quire.PositExpRange}.modularSub(exp) ;
	auto frac = lzocshifted.template slice<
			allsize - 2,
			allsize - (PositDim<N, WES>::WF + 1)
		>();
	auto guard = lzocshifted.template get<allsize - PositDim<N, WES>::WF - 2>();
	auto stickycomp = lzocshifted.template slice<
			allsize - PositDim<N, WES>::WF - 3,
			0
		>();
	auto sticky = stickycomp.or_reduction() | lower_sticky;

	auto isSpecial = isZero | fin_overflow | underflow | quire.getIsNaR();

	auto fin_sticky = sticky & isSpecial.invert();
	auto fin_guard = guard & isSpecial.invert();
	/*
	cerr << to_string(fin_overflow) << endl;
	cerr << to_string(isZero) << endl;
	cerr << to_string(underflow) << endl;
	*/

	constexpr unsigned int overflow_neg = 2*PositDim<N, WES>::EXP_BIAS;
	constexpr unsigned int overflow_pos = overflow_neg - 1;

	Wrapper<logSize, false> overflow_neg_wrap{{ overflow_neg }};
	Wrapper<logSize, false> overflow_pos_wrap{{ overflow_pos }};

	auto overflow_exp_val = Wrapper<logSize, false>::mux(sign, overflow_neg_wrap, overflow_pos_wrap);

	auto fin_exp = Wrapper<logSize, false>::mux(
				high_overflow,
				overflow_exp_val,
				Wrapper<logSize, false>::mux(
					isZero,
					{0},
					Wrapper<logSize, false>::mux(
						underflow,
						sign.invert().template leftpad<logSize>(),
						biased_exp
					)
				)
		);

	auto implicit_bit = sign.invert() & isZero.invert();
	auto fin_frac = Wrapper<PositDim<N, WES>::WF, false>::mux(isSpecial, {0}, frac);
	/*
	cerr << to_string(fin_guard) << endl;
	cerr << to_string(fin_sticky) << endl;
	cerr << to_string(fin_exp) << endl;
	cerr << to_string(sign) << endl;
	cerr << to_string(implicit_bit) << endl;
	cerr << to_string(fin_frac) << endl;
	*/
	return PositIntermediateFormat<N, WES, Wrapper, false>(
			fin_guard,
			fin_sticky,
			quire.getIsNaR(),
			fin_exp,
			sign,
			implicit_bit,
			fin_frac
		);
}

