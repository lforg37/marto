#pragma once
#include "posit_dim.hpp"
#include "lzoc_shifter.hpp"


template<int N> 
ap_uint<PositDim<N>::WE> getExponent(
		ap_uint<get2Power(N)> range_count,
		ap_uint<N-3> shifted_fraction,
		ap_uint<1> sign,
		typename enable_if<PositDim<N>::HAS_ES>::type* dummy = 0
	) 
{
	ap_uint<PositDim<N>::WES> es = shifted_fraction.range(N-4, N - 3 - PositDim<N>::WES); 
	ap_int<PositDim<N>::WES> signed_ext_sign = (ap_int<1>) sign;
	ap_int<PositDim<N>::WES> ext_sign = signed_ext_sign;

	ap_uint<PositDim<N>::WES> decoded_es = es^ext_sign;

	return range_count.concat(decoded_es);
}

template<int N> 
ap_uint<PositDim<N>::WE> getExponent(
		ap_uint<get2Power(N)> range_count,
		ap_uint<N-3> shifted_fraction,
		ap_uint<1> sign,
		typename enable_if<not PositDim<N>::HAS_ES>::type* dummy = 0
	) 
{
	return range_count;
}

template<int N>
PositValue<N> posit_decoder(PositEncoding<N> positN)
{
	//Sign bit
	ap_uint<1> s = positN[N-1];
	//First regime bit
	ap_uint<1> count_type = positN[N-2];
	//Remainder
	ap_uint<N-2> remainder = positN.range(N-3, 0);

	//Add !count_type to the end of the input
	ap_int<2> spad = (ap_int<1>) (not count_type);
	ap_uint<2> pad = spad;
	ap_uint<N> input_shift = remainder.concat(pad);

	ap_uint<1> ZeroNAR = remainder.or_reduce();
	ap_uint<1> isNAR = ZeroNAR and s;
	
	ap_uint<1> implicit_bit = not(s) and not(ZeroNAR);

	constexpr int logN = get2Power(N);
	auto lzoc_shifted = lzoc_shifter<logN>(input_shift, count_type);
	ap_uint<logN> RangeCount = lzoc_shifted.range(N+logN-1, N);
	ap_uint<N-3> usefulBits = lzoc_shifted(N-2, 2);
	ap_uint<N-3-PositDim<N>::WES> fraction = usefulBits(N-4-PositDim<N>::WES, 0);
	
	ap_int<1> neg_count = not(s^count_type);
	ap_int<logN> signed_extended_neg_count = neg_count;
	ap_uint<logN> extended_neg_count = signed_extended_neg_count;

	auto exponent = getExponent<N>(RangeCount, usefulBits, s);	
	return PositValue<N>(
			isNAR, 
			exponent,
			s,
			implicit_bit,
			fraction
		);
}
