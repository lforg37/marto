#ifndef QUIRE_TO_POSIT_HPP
#define QUIRE_TO_POSIT_HPP

#include "lzoc_shifter.hpp"
#include "posit_dim.hpp"

template<int N>
PositValue<N> quire_to_posit(Quire<N> quire)
{
	constexpr int size_range = 2*Quire<N>::PositRangeOffset;
	constexpr int logSize = ceilLog2(size_range);
	constexpr int allsize = ceil2Power(size_range);
	constexpr int padd_width = allsize - size_range;

	ap_int<1> sign = quire.getSignBit();

	ap_uint<padd_width> upper_low_bits = quire.range(
			Quire<N>::PositRangeOffset - 1, 
			Quire<N>::PositRangeOffset - padd_width
		);

	ap_uint<Quire<N>::PositRangeOffset - padd_width> lower_low_bits = quire.range(
			Quire<N>::PositRangeOffset - padd_width - 1,
			0
		);

	ap_uint<1> lower_sticky = not(lower_low_bits.or_reduce());
	ap_uint<1> upper_low_null = not(upper_low_bits.or_reduce());

	ap_uint<1> low_bit_is_null = lower_sticky and upper_low_null;

	ap_uint<2*Quire<N>::PositRangeOffset> middle_bits = quire.range(
			Quire<N>::PositRangeOffset,
			3 * Quire<N>::PositRangeOffset - 1
		);

	ap_int<2*Quire<N>::PositRangeOffset> middle_s_ext = not sign;
	ap_uint<2*Quire<N>::PositRangeOffset> underflow_base = middle_s_ext xor middle_bits;
	ap_uint<1> middle_void_flag = underflow_base.or_reduce();
	ap_uint<1> middle_is_null = not(middle_bits.or_reduce());

	ap_uint<1> uperincludedbound = middle_bits[2*Quire<N>::PositRangeOffset - 1];

	constexpr int remainingsize = 
		PositDim<N>::ExtQuireSize - 2 - 4 * Quire<N>::PositRangeOffset;

	ap_uint<remainingsize> uppercarry = quire.range(PositDim<N>::ExtQuireSize - 3,
			4*Quire<N>::PositRangeOffset
		); 

	ap_int<remainingsize> upper_ext_sign = sign;
	ap_uint<remainingsize> overflow = upper_ext_sign xor uppercarry;
	ap_uint<1> high_overflow = overflow.or_reduce();

	ap_uint<1> fin_overflow = 
		high_overflow or (sign and middle_is_null and low_bit_is_null);

	ap_uint<1> isZero = (not sign) and 
						(not fin_overflow) and 
						middle_is_null and 
						low_bit_is_null;

	ap_uint<1> underflow = 	not(underflow)
		and not(overflow)
		and middle_void_flag 
		and not(sign and low_bit_is_null);

	
	ap_uint<allsize> padded_mid_bits = middle_bits.concat(
			( (ap_uint<padd_width>)
			quire.range(
					Quire<N>::PositRangeOffset - 1, 
					Quire<N>::PositRangeOffset - padd_width
				)
			)
			);

	auto lzocshifted = lzoc_shifter<logSize>(padded_mid_bits, sign);
	ap_uint<logSize> exp = lzocshifted.range(logSize + allsize - 1, allsize);

	ap_uint<logSize> biased_exp = ap_uint<logSize>{size_range} - exp;
	ap_uint<PositDim<N>::WF> frac = lzocshifted.range(allsize - 2, allsize - (PositDim<N>::WF + 1));
	ap_uint<1> guard = lzocshifted[allsize - PositDim<N>::WF - 2];
	ap_uint<allsize -( PositDim<N>::WF + 2)> stickycomp = lzocshifted.range(
			allsize - PositDim<N>::WF - 3, 
			0
		);
	ap_uint<1> sticky = ((ap_uint<1>) stickycomp.or_reduce()) or lower_sticky;

	ap_uint<1> isSpecial = isZero or fin_overflow or underflow or quire.getIsNaR();

	ap_uint<1> fin_sticky = sticky and (not isSpecial);
	ap_uint<1> fin_guard = guard and (not isSpecial);

	ap_uint<logSize> fin_exp;
	if (overflow) {
		fin_exp = 2*PositDim<N>::EXP_BIAS - 1 - sign;
	} else if (isZero) {
		fin_exp = 0;
	} else if (underflow)
	{
		fin_exp = not sign;
	} else {
		fin_exp = biased_exp;
	}

	ap_uint<1> implicit_bit = (not sign) and (not isZero); 
	ap_uint<PositDim<N>::WF> fin_frac;
	if (isSpecial) {
		fin_frac = 0;
	} else {
		fin_frac = frac;
	}
	
	return PositValue<N>(
			fin_guard,
			fin_sticky,
			quire.getIsNaR(), 
			fin_exp,
			sign,
			implicit_bit,
			fin_frac
		);
}
#endif
