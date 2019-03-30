#ifndef VALUE_PROD_CONVERSIONS_HPP
#define VALUE_PROD_CONVERSIONS_HPP

#include <cstdio> 

#include "posit_dim.hpp"

template<int N, int WES>
PositProd<N, WES> PositValue_to_PositProd(PositValue<N, WES> val)
{
	#pragma HLS INLINE
    ap_uint<PositValue<N, WES>::FractionSize + 1> signed_val_frac = val.getSignificand();

    ap_uint<PositProd<N, WES>::SignificandSize> significand = signed_val_frac.concat(
            ap_uint<PositValue<N, WES>::FractionSize+1>(0)
		);


    auto expt = val.getExp();

    ap_uint<PositProd<N, WES>::ExpSize> exponent;
   	if (val.isZero()) {
		exponent = 0;	
	} else {
        exponent = static_cast<ap_uint<PositProd<N, WES>::ExpSize> >(val.getExp()) +
            ap_uint<PositProd<N, WES>::ExpSize>{PositDim<N, WES>::EXP_BIAS - 1};
	}

    return PositProd<N, WES>(val.getIsNaR(), exponent, val.getSignBit(), significand);
}

template<int N, int WES>
PositValue<N, WES> PositProd_to_PositValue(PositProd<N, WES> val)
{
	#pragma HLS INLINE
	ap_uint<1> isMinPos, isMaxPos;
	ap_uint<1> isNaR = val.getIsNaR();
    ap_uint<PositProd<N, WES>::SignificandSize> fraction = val.getSignificand();
    ap_uint<1> implicitBit = fraction[PositProd<N, WES>::SignificandSize-1];
    ap_uint<PositValue<N, WES>::FractionSize> resultFraction = fraction.range(PositProd<N, WES>::SignificandSize-1-1, PositProd<N, WES>::SignificandSize-1-1 -(PositValue<N, WES>::FractionSize)+1);
    ap_uint<1> resultGuardBit = fraction[PositProd<N, WES>::SignificandSize-1-1 -(PositValue<N, WES>::FractionSize)+1-1];
    ap_uint<1> resultStickyBit = not(fraction.range(PositProd<N, WES>::SignificandSize-1-1 -(PositValue<N, WES>::FractionSize)+1-1-1,0) == 0);

	// ap_uint<1> isZero = not(((ap_uint<4>) val.getSignificand().range(PositDim<N>::ProdSignificandSize - 1, PositDim<N>::ProdSignificandSize - 4)).or_reduce());
	ap_uint<1> isZero = not(val.getSignBit()) and not(implicitBit) and val.getExp()==0;

    ap_int<PositProd<N, WES>::ExpSize+1> exp;
	
    ap_uint<PositProd<N, WES>::ExpSize> expt = val.getExp();
	// printApUint(expt);	
	// printApUint(bias);	
	// printApUint(isZero);	

	if (isZero) {
		exp = 0;
	} else {
        exp	= static_cast<ap_int<PositProd<N, WES>::ExpSize+1> >(val.getExp()-(PositDim<N, WES>::EXP_BIAS - 1));
	}
	// printApInt(exp);	


	if(exp<0 and not(isNaR)){
		isMinPos = 1;
	} else {
		isMinPos = 0;
	}
    if(exp >= (2*PositDim<N, WES>::EXP_BIAS) and not(isNaR)){
		isMaxPos = 1;
	}
	else{
		isMaxPos = 0;
	}

	// fprintf(stderr, "isminpos: ");
	// printApUint(isMinPos);
	// fprintf(stderr, "ismaxpos: ");
	// printApUint(isMaxPos);


	if(isMinPos and val.getSignBit() ){
        return PositValue<N, WES>::getMinNeg();
	}
	else if(isMinPos and not(val.getSignBit())){
        return PositValue<N, WES>::getMinPos();
	}
	else if(isMaxPos and val.getSignBit()){
        return PositValue<N, WES>::getMaxNeg();
	}
	else if(isMaxPos and not(val.getSignBit())){
        return PositValue<N, WES>::getMaxPos();
	}
	else {
        return PositValue<N, WES>(
				resultGuardBit,
				resultStickyBit,
				isNaR,
                exp.range(PositValue<N, WES>::ExpSize-1,0), //Warning : biased exp
				val.getSignBit(),
				implicitBit,
				resultFraction
            );
	}
}
#endif
