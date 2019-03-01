#ifndef POSIT_MUL_HPP
#define POSIT_MUL_HPP

#include "posit_dim.hpp"

using namespace std;

template <int N> PositProd<N> posit_mul(PositValue<N> in1, PositValue<N> in2)
{
	#pragma HLS INLINE
    ap_uint<1> isNar = in1.getIsNaR() | in2.getIsNaR();
    ap_int<1> isZero = in1.isZero() or in2.isZero();

    // Compute the significand
    ap_uint<PositDim<N>::ProdSignificandSize + 2> significand = 
		in1.getSignedSignificand() * in2.getSignedSignificand();

	ap_uint<1> sign = significand[PositDim<N>::ProdSignificandSize + 1];
	ap_uint<1> exbit = significand[PositDim<N>::ProdSignificandSize - 1];
	ap_uint<1> neg_neg_2power = sign xor significand[PositDim<N>::ProdSignificandSize];

	ap_uint<1> needs_shift = (exbit xor sign) or neg_neg_2power;

	ap_uint<PositDim<N>::ProdSignificandSize> fin_significand;
	if (needs_shift) {
		ap_uint<1> first_bit = (significand[PositDim<N>::ProdSignificandSize - 1]) 
			or neg_neg_2power;

		ap_uint<PositDim<N>::ProdSignificandSize - 1> last_bits = 
			significand.range(PositDim<N>::ProdSignificandSize - 2, 0);

		fin_significand = first_bit.concat(last_bits);
	} else {
		fin_significand = ((ap_uint<PositDim<N>::ProdSignificandSize - 1>) 
				significand.range(PositDim<N>::ProdSignificandSize - 2, 0)).concat(ap_uint<1>{0}); 
	}

    // Compute the exponent
    ap_uint<PositDim<N>::ProdExpSize> exponent;
    if (isZero) {
		exponent = 0;
    } else {
		exponent = in1.getExp() + in2.getExp() - not(needs_shift) + neg_neg_2power;
    }

    return PositProd<N>(
		isNar, 
		exponent, 
		sign,
		fin_significand
	);
}
#endif
