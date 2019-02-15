#ifndef VALUE_PROD_CONVERSIONS_HPP
#define VALUE_PROD_CONVERSIONS_HPP

#include "posit_dim.hpp"

template<int N>
PositProd<N> PositValue_to_PositProd(PositValue<N> val) 
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

template<int N>
PositValue<N> PositProd_to_PositValue(PositProd<N> val) 
{

	ap_int<PositDim<N+1>::ProdExpSize> exp = (ap_int<PositDim<N+1>::ProdExpSize>)val.getExp()-PositDim<N>::EXP_BIAS;
	ap_uint<1> isMinPos, isMaxPos;
	ap_uint<1> isNaR = val.getIsNaR();

	if(exp<0 and not(isNaR)){
		isMinPos = 1;
		isMaxPos = 0;
	}
	else if(exp >= (2*PositDim<N>::EXP_BIAS) and not(isNaR)){
		isMinPos = 0;
		isMaxPos = 1;
	}
	else{
		isMinPos = 0;
		isMaxPos = 0;
	}
	// printApInt(exp);

	ap_uint<PositDim<N>::ProdSignificandSize> fraction = val.getSignificand();
	// printApUint(fraction);

	if(isMinPos and val.getSignBit() ){
		return PositValue<N>::getMinNeg();
	}
	else if(isMinPos and not(val.getSignBit())){
		return PositValue<N>::getMinPos();
	}
	else if(isMaxPos and val.getSignBit()){
		return PositValue<N>::getMaxNeg();
	}
	else if(isMaxPos and not(val.getSignBit())){
		return PositValue<N>::getMaxPos();
	}
	else{
		return PositValue<N>(
				// ap_uint<1> guard,
				0,
				// ap_uint<1> sticky,
				0,
				isNaR,
				// ap_uint<PositDim<N>::WE> exp, //Warning : biased exp
				exp.range(PositDim<N>::WE-1,0), //Warning : biased exp
				val.getSignBit(),
				val.getSignificand()[PositDim<N>::ProdSignificandSize-4],
				// ap_uint<PositDim<N>::WF> fraction);
				val.getSignificand().range(PositDim<N>::ProdSignificandSize-5, PositDim<N>::ProdSignificandSize-5 -(PositDim<N>::WF)+1)
				);
	}
}


#endif
