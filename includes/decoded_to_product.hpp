#ifndef DECODED_TO_PRODUCT_HPP
#define DECODED_TO_PRODUCT_HPP

#include "posit_dim.hpp"

template<int N>
PositProd<N> decoded_to_product(PositValue<N> val) 
{
	ap_int<PositDim<N>::WF + 4> signed_val_frac = val.getSignedSignificand();
	ap_uint<PositDim<N>::WF + 4> sign_extended_frac = signed_val_frac;

	ap_uint<PositDim<N>::ProdSignificandSize> significand = sign_extended_frac.concat(
			ap_uint<PositDim<N>::WF>(0)
		);

	ap_uint<PositDim<N>::ProdExpSize> exponent = 
		((ap_uint<PositDim<N>::ProdExpSize>) val.getExp()) + 
		ap_uint<PositDim<N>::ProdExpSize>(PositDim<N>::EXP_BIAS);

	return PositProd<N>(val.getIsNaR(), exponent, significand);
}

#endif
