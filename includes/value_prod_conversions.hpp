#ifndef VALUE_PROD_CONVERSIONS_HPP
#define VALUE_PROD_CONVERSIONS_HPP

#include <cstdio> 

#include "posit_dim.hpp"

template<int N>
PositProd<N> PositValue_to_PositProd(PositValue<N> val) 
{
	#pragma HLS INLINE
	ap_uint<PositDim<N>::WF + 1> signed_val_frac = val.getSignificand();

	ap_uint<PositDim<N>::ProdSignificandSize> significand = signed_val_frac.concat(
			ap_uint<PositDim<N>::WF+1>(0)
		);


	ap_uint<PositDim<N>::WE> expt = val.getExp();

	ap_uint<PositDim<N>::ProdExpSize> exponent;
   	if (val.isZero()) {
		exponent = 0;	
	} else {
		exponent = ((ap_uint<PositDim<N>::ProdExpSize>) val.getExp()) + 
		ap_uint<PositDim<N>::ProdExpSize>(PositDim<N>::EXP_BIAS - 1);
	}

	return PositProd<N>(val.getIsNaR(), exponent, val.getSignBit(), significand);
}

extern template PositProd<8> PositValue_to_PositProd(PositValue<8> val); 
extern template PositProd<16> PositValue_to_PositProd(PositValue<16> val); 
extern template PositProd<32> PositValue_to_PositProd(PositValue<32> val); 
extern template PositProd<64> PositValue_to_PositProd(PositValue<64> val); 

template<int N>
PositValue<N> PositProd_to_PositValue(PositProd<N> val) 
{
	#pragma HLS INLINE
	ap_uint<1> isMinPos, isMaxPos;
	ap_uint<1> isNaR = val.getIsNaR();
	ap_uint<PositDim<N>::ProdSignificandSize> fraction = val.getSignificand();
	ap_uint<1> implicitBit = fraction[PositDim<N>::ProdSignificandSize-1];
	ap_uint<PositDim<N>::WF> resultFraction = fraction.range(PositDim<N>::ProdSignificandSize-1-1, PositDim<N>::ProdSignificandSize-1-1 -(PositDim<N>::WF)+1);
	ap_uint<1> resultGuardBit = fraction[PositDim<N>::ProdSignificandSize-1-1 -(PositDim<N>::WF)+1-1];
	ap_uint<1> resultStickyBit = not(fraction.range(PositDim<N>::ProdSignificandSize-1-1 -(PositDim<N>::WF)+1-1-1,0) == 0);

	// ap_uint<1> isZero = not(((ap_uint<4>) val.getSignificand().range(PositDim<N>::ProdSignificandSize - 1, PositDim<N>::ProdSignificandSize - 4)).or_reduce());
	ap_uint<1> isZero = not(val.getSignBit()) and not(implicitBit) and val.getExp()==0;

	ap_int<PositDim<N>::ProdExpSize+1> exp;
	
	ap_uint<PositDim<N>::ProdExpSize> expt = val.getExp();
	// printApUint(expt);	
	// printApUint(bias);	
	// printApUint(isZero);	

	if (isZero) {
		exp = 0;
	} else {
		exp	= (ap_int<PositDim<N>::ProdExpSize+1>)val.getExp()-(PositDim<N>::EXP_BIAS - 1);
	}
	// printApInt(exp);	


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

	// fprintf(stderr, "isminpos: ");
	// printApUint(isMinPos);
	// fprintf(stderr, "ismaxpos: ");
	// printApUint(isMaxPos);


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
				resultGuardBit,
				resultStickyBit,
				isNaR,
				exp.range(PositDim<N>::WE-1,0), //Warning : biased exp
				val.getSignBit(),
				implicitBit,
				resultFraction
				);
	}
}

extern template PositValue<8> PositProd_to_PositValue(PositProd<8> val);
extern template PositValue<16> PositProd_to_PositValue(PositProd<16> val);
extern template PositValue<32> PositProd_to_PositValue(PositProd<32> val);
extern template PositValue<64> PositProd_to_PositValue(PositProd<64> val);
#endif
