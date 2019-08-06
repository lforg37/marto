#pragma once
#include <iostream>

#include "posit_dim.hpp"

#define S_WF PositIntermediateFormat<N, WES>::FractionSize
#define S_WE PositIntermediateFormat<N, WES>::ExpSize
#define S_WES WES
#define K_SIZE (S_WE-S_WES)

template<int N, int WES>
PositEncoding<N, WES> posit_encoder(PositIntermediateFormat<N, WES> positValue)
{
    ap_uint<S_WE> expWoBias = positValue.getExp() - PositDim<N, WES>::EXP_BIAS;

	ap_uint<1> sign = positValue.getSignBit();
	ap_uint<S_WES> es = expWoBias.range(S_WES-1,0) ^ sign;
	ap_int<K_SIZE> k = expWoBias.range(S_WE-1,S_WES);
	ap_uint<S_WF> exactSignificand = positValue.getFraction();
	ap_uint<N-1-2-S_WES> significand = exactSignificand.range(S_WF-1, S_WF-1-(N-1-2-S_WES)+1);

	ap_uint<S_WES+N-1-2-S_WES> esAndSignificand = es.concat(significand);
    ap_uint<2> zero_one{1};
    ap_uint<2> one_zero{2};

	ap_int<2+S_WES+N-1-2-S_WES> reverseBitAndEsAndSignificand;
	if((k[K_SIZE-1] == 1)^sign){
		reverseBitAndEsAndSignificand = zero_one.concat(esAndSignificand);
	}
	else{
		reverseBitAndEsAndSignificand = one_zero.concat(esAndSignificand);
	}

	ap_uint<K_SIZE-1> absK;

	if(k[K_SIZE-1] == 1){
		absK = ~k;
	}
	else{
		absK = k;
	}

    ap_int<N> reverseBitAndEsAndSignificandAndGuardBit = reverseBitAndEsAndSignificand.concat(positValue.getGuardBit());
	// printApInt(readyToShift);

    ap_uint<N+1> shifted = shifter_sticky<N, K_SIZE-1, true>(
                reverseBitAndEsAndSignificandAndGuardBit,
                absK,
                reverseBitAndEsAndSignificandAndGuardBit[N-1]
        ); //TODO rajouter le fillbit

	// printApUint(shifted);

    ap_uint<N-1> unroundedResult =  shifted.range(N, 2);

	// printApUint(unroundedResult);

    ap_uint<1> guard = shifted[1];
    ap_uint<1> sticky = shifted[0] or positValue.getStickyBit();

	// printApUint(guard);
	// printApUint(sticky);


	ap_uint<1> roundingBit = (guard and not(sticky) and unroundedResult[0]) or (guard and sticky);

    ap_uint<N-1> roundedResult = unroundedResult + roundingBit;

	ap_uint<N> normalOutput = sign.concat(roundedResult);
	ap_uint<N-1> zero = 0;
	ap_uint<1> isNaRBit = positValue.getIsNaR();
	ap_uint<N> specialCasesValue = isNaRBit.concat(zero);
	return (((not positValue.getSignBit()) and (not positValue.getImplicitBit())) or isNaRBit) ? specialCasesValue : normalOutput;
}
