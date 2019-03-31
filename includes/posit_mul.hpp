#pragma once

#include "posit_dim.hpp"

using namespace std;

template <int N, int WES> PositProd<N, WES> posit_mul(PositValue<N, WES> in1, PositValue<N, WES> in2)
{
	#pragma HLS INLINE
    ap_uint<1> isNar = in1.getIsNaR() | in2.getIsNaR();
    ap_int<1> isZero = in1.isZero() or in2.isZero();

    // Compute the significand
    ap_uint<PositProd<N, WES>::SignificandSize + 2> significand =
		in1.getSignedSignificand() * in2.getSignedSignificand();

    ap_uint<1> sign = significand[PositProd<N, WES>::SignificandSize + 1];
    ap_uint<1> exbit = significand[PositProd<N, WES>::SignificandSize - 1];
    ap_uint<1> neg_neg_2power = sign xor significand[PositProd<N, WES>::SignificandSize];

	ap_uint<1> needs_shift = (exbit xor sign) or neg_neg_2power;

    ap_uint<PositProd<N, WES>::SignificandSize> fin_significand;
    ap_uint<PositProd<N, WES>::SignificandSize - 1> last_bits =
        significand.range(PositProd<N, WES>::SignificandSize - 2, 0);
	if (needs_shift) {
        ap_uint<1> first_bit = (significand[PositProd<N, WES>::SignificandSize - 1])
			or neg_neg_2power;
		fin_significand = first_bit.concat(last_bits);
	} else {
		fin_significand = last_bits.concat(ap_uint<1>{0}); 
	}

    // Compute the exponent
    ap_uint<PositProd<N, WES>::ExpSize> exponent;
    if (isZero) {
		exponent = 0;
    } else {
		exponent = in1.getExp() + in2.getExp() - not(needs_shift) + neg_neg_2power;
    }

    return PositProd<N, WES>(
		isNar, 
		exponent, 
		sign,
		fin_significand
	);
}
