#ifndef POSIT_MUL_HPP
#define POSIT_MUL_HPP

#include "posit_dim.hpp"

using namespace std;

template<int N>
PositProd<N> posit_mul(PositValue<N> in1, PositValue<N> in2)
{
	ap_uint<1> isNar = in1.getIsNaR() | in2.getIsNaR();
	ap_int<1> isZero = in1.isZero() or in2.isZero();

	//Compute the significand
	ap_uint<PositDim<N>::ProdSignificandSize> significand = 
		in1.getSignedSignificand() * in2.getSignedSignificand();

	//Compute the exponent
	ap_uint<PositDim<N>::ProdExpSize> exponent = 
		in1.getExp() + in2.getExp();

	ap_int<PositDim<N>::ProdExpSize> zeroMask = not isZero;
	ap_uint<PositDim<N>::ProdExpSize> fin_exp = exponent and zeroMask;

	return PositProd<N>(isNar, fin_exp, significand);
}
#endif
