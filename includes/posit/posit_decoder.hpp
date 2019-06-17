#pragma once
#include "posit_dim.hpp"
#include "marto/bitvector.hpp"

template<int N, int WES>
ap_uint<PositDim<N, WES>::WE> getExponent(
		ap_uint<Static_Val<N>::_log2 + 1> range_count,
		ap_uint<N-3> shifted_fraction,
		ap_uint<1> sign,
        typename enable_if<PositDim<N, WES>::HAS_ES>::type* = 0
	) 
{
	#pragma HLS INLINE
    ap_uint<WES> es = shifted_fraction.range(N-4, N - 3 - WES);
    ap_int<WES> signed_ext_sign = static_cast<ap_int<1> >(sign);
    ap_int<WES> ext_sign = signed_ext_sign;

    ap_uint<WES> decoded_es = es^ext_sign;

	return range_count.concat(decoded_es);
}

template<int N, int WES>
ap_uint<PositDim<N, WES>::WE> getExponent(
		ap_uint<Static_Val<N>::_log2 + 1> range_count,
        ap_uint<N-3>,
        ap_uint<1>,
        typename enable_if<not PositDim<N, WES>::HAS_ES>::type* = 0
	) 
{
	#pragma HLS INLINE
	return range_count;
}

template<int N, int WES>
PositIntermediateFormat<N, WES> posit_decoder(PositEncoding<N, WES> positN)
{
	#pragma HLS INLINE
	//Sign bit
	ap_uint<1> s = positN[N-1];
	//First regime bit
	ap_uint<1> count_type = positN[N-2];
	//Remainder
	ap_uint<N-2> remainder = positN.range(N-3, 0);

	//Add 0 to the end of the input
	ap_uint<2> pad = 0;
	ap_uint<N> input_shift = remainder.concat(pad);

	ap_uint<1> zero_NAR = not(remainder.or_reduce() or count_type);
	ap_uint<1> is_NAR = zero_NAR and s;
	ap_uint<1> is_zero = zero_NAR and not s;
	
	ap_uint<1> implicit_bit = not(s) and not(zero_NAR);

	constexpr int logN = Static_Val<N>::_log2;
	auto lzoc_shifted = lzoc_shifter<logN>(input_shift, count_type);
	ap_uint<logN + 1> RangeCount = lzoc_shifted.range(N+logN-1, N);
	ap_uint<N-3> usefulBits = lzoc_shifted.range(N-2, 2);
    ap_uint<N-3-WES> fraction = usefulBits(N-4-WES, 0);
	
	ap_int<1> neg_count = not(s^count_type);
	ap_int<logN + 1> signed_extended_neg_count = neg_count;
	ap_uint<logN + 1> extended_neg_count = signed_extended_neg_count;
	ap_uint<logN + 1> comp2_range_count = RangeCount^extended_neg_count;

    ap_uint<PositIntermediateFormat<N, WES>::ExpSize> exponent = getExponent<N, WES>(
			comp2_range_count,
		   	usefulBits,
		   	s
		);
   ap_uint<PositIntermediateFormat<N, WES>::ExpSize> biased_exp = exponent + ap_uint<PositIntermediateFormat<N, WES>::ExpSize>(PositDim<N, WES>::EXP_BIAS);

   ap_int<1> is_not_zero = not is_zero;
   ap_int<PositIntermediateFormat<N, WES>::ExpSize> extended_is_not_zero = is_not_zero;
   ap_uint<PositIntermediateFormat<N, WES>::ExpSize> unsigned_ext_is_not_zero = extended_is_not_zero;
   ap_uint<PositIntermediateFormat<N, WES>::ExpSize> final_biased_exp = biased_exp & unsigned_ext_is_not_zero;

    return PositIntermediateFormat<N, WES>(
			is_NAR, 
			final_biased_exp,
			s,
			implicit_bit,
			fraction
		);
}