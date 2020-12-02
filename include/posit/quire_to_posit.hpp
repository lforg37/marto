#pragma once

#include <iostream>

#include "posit_dim.hpp"
#include "primitives/lzoc_shifter.hpp"

#ifdef POSIT_QUIRETOPIF_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif

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
	auto exp = Wrapper<logSize, false>{{PositDim<N, WES>::EMax}}.modularSub(lzocshifted.lzoc);

	auto frac = lzocshifted.shifted.template slice<
			allsize - 2,
			allsize - (PositDim<N, WES>::WF + 1)
		>();
	auto guard = lzocshifted.shifted.template get<allsize - PositDim<N, WES>::WF - 2>();
	auto stickycomp = lzocshifted.shifted.template slice<
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

	constexpr unsigned int overflow_neg = PositDim<N, WES>::EMax-1;
	Wrapper<logSize, false> overflow_neg_wrap{{ overflow_neg }};

	constexpr unsigned int underflow_neg = ((overflow_neg + 1)^ (((1) << PositDim<N, WES>::WE)-1)) ;
	Wrapper<logSize, false> underflow_neg_wrap{{underflow_neg}};


	auto overflow_base = Wrapper<logSize, false>::mux(high_overflow, overflow_neg_wrap, underflow_neg_wrap);
	auto over_flow_exp = overflow_base.modularAdd(sign.invert().template leftpad<logSize>());

	auto overflow_fin_exp = Wrapper<logSize, false>::mux(
				high_overflow | underflow,
				over_flow_exp,
				exp
		);

	auto isZeroMask = Wrapper<logSize, false>::generateSequence(isZero.invert());
	auto fin_exp = overflow_fin_exp & isZeroMask;


	auto implicit_bit = sign.invert() & isZero.invert();
	auto isSpecialMask = Wrapper<PositDim<N, WES>::WF, false>::generateSequence(isSpecial.invert());
	auto fin_frac = frac & isSpecialMask;

#ifdef POSIT_QUIRETOPIF_DEBUG
	cerr << "=== QUIRE to PIF ===" << endl;
	cerr << "sign: " << to_string(sign) << endl;
	cerr << "upper_low_bits: " << to_string(upper_low_bits) << endl;
	cerr << "lower_low_bits: " << to_string(lower_low_bits) << endl;
	cerr << "lower_sticky: " << to_string(lower_sticky) << endl;
	cerr << "upper_low_null: " << to_string(upper_low_null) << endl;
	cerr << "low_bit_is_null: " << to_string(low_bit_is_null) << endl;
	cerr << "middle_bits: " << to_string(middle_bits) << endl;
	cerr << "middle_s_ext: " << to_string(middle_s_ext) << endl;
	cerr << "underflow_base: " << to_string(underflow_base) << endl;
	cerr << "middle_void_flag: " << to_string(middle_void_flag) << endl;
	cerr << "middle_is_null: " << to_string(middle_is_null) << endl;
	cerr << "uppercarry: " << to_string(uppercarry) << endl;
	cerr << "upper_ext_sign: " << to_string(upper_ext_sign) << endl;
	cerr << "overflow: " << to_string(overflow) << endl;
	cerr << "high_overflow: " << to_string(high_overflow) << endl;
	cerr << "fin_overflow: " << to_string(fin_overflow) << endl;
	cerr << "isZero: " << to_string(isZero) << endl;
	cerr << "underflow: " << to_string(underflow) << endl;
	cerr << "lzocshifted: " << to_string(lzocshifted) << endl;
	cerr << "exp: " << to_string(exp) << endl;
	cerr << "frac: " << to_string(frac) << endl;
	cerr << "guard: " << to_string(guard) << endl;
	cerr << "stickycomp: " << to_string(stickycomp) << endl;
	cerr << "sticky: " << to_string(sticky) << endl;
	cerr << "isSpecial: " << to_string(isSpecial) << endl;
	cerr << "fin_sticky: " << to_string(fin_sticky) << endl;
	cerr << "fin_guard: " << to_string(fin_guard) << endl;
	cerr << "overflow_base: " << to_string(overflow_base) << endl;
	cerr << "over_flow_exp: " << to_string(over_flow_exp) << endl;
	cerr << "overflow_fin_exp: " << to_string(overflow_fin_exp) << endl;
	cerr << "isZeroMask: " << to_string(isZeroMask) << endl;
	cerr << "fin_exp: " << to_string(fin_exp) << endl;
	cerr << "implicit_bit: " << to_string(implicit_bit) << endl;
	cerr << "isSpecialMask: " << to_string(isSpecialMask) << endl;
	cerr << "fin_frac: " << to_string(fin_frac) << endl;
	cerr << "===============================" << endl;
#endif

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

