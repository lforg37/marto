#ifndef POSIT_MUL_HPP
#define POSIT_MUL_HPP

#include "posit_dim.hpp"

using namespace std;

template<int N>
PositProd<N> posit_mul(PositValue<N> in1, PositValue<N> in2)
{
	ap_uint<1> isNar = in1.getIsNaR() | in2.getIsNaR();

	//Compute the significand
	ap_uint<PositDim<N>::ProdSignificandSize> significand = 
		in1.getSignedSignificand() * in2.getSignedSignificand();

	//Compute the exponent
	ap_uint<PositDim<N>::ProdExpSize> exponent = 
		in1.getExp() + in2.getExp();

	return PositProd<N>(isNar, exponent, significand);
}
#endif
