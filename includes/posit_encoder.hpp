#pragma once
#include "posit_dim.hpp"

#define S_WF PositDim<N>::WF
#define S_WE PositDim<N>::WE
#define S_WES PositDim<N>::WES
#define K_SIZE (S_WE-S_WES)

template<int N>
PositEncoding<N> posit_encoder(PositValue<N> positValue)
{
	// ap_uint<1> s_in = fp[WE+WF];
	// ap_uint<WE> fp_exp_in = fp.range(WE+WF-1, WF);
	// ap_uint<WF> fp_mant_in = fp.range(WF-1, 0);

	// fp_exp_in -= FP_BIAS;

	ap_uint<S_WE> expWoBias = positValue.getExp() - PositDim<N>::EXP_BIAS;

	ap_uint<1> sign = positValue.getSignBit();
	ap_uint<S_WES> es = expWoBias.range(S_WES-1,0) ^ sign;
	ap_int<K_SIZE> k = expWoBias.range(S_WE-1,S_WES);
	ap_uint<S_WF> exactSignificand = positValue.getSignificandWoImp();
	ap_uint<N-1-2-S_WES> significand = exactSignificand.range(S_WF-1, S_WF-1-(N-1-2-S_WES)+1);

	ap_uint<S_WES+N-1-2-S_WES> esAndSignificand = es.concat(significand);
	ap_uint<2> zero_one = 0b01;
	ap_uint<2> one_zero = 0b10;

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

	ap_uint<2+S_WES+N-1-2-S_WES> shiftedResverseBitAndEsAndSignificand = reverseBitAndEsAndSignificand >> absK;

	// ap_uint<2+ES+N-1-2-ES>	shiftedResverseBitAndEsAndSignificandComp;
	// if(s_in == 1){
	// 	shiftedResverseBitAndEsAndSignificandComp = ~shiftedResverseBitAndEsAndSignificand+1;
	// }
	// else{
	// 	shiftedResverseBitAndEsAndSignificandComp = shiftedResverseBitAndEsAndSignificand;	
	// }

	ap_uint<N-1> trunctatedEsAndSignificand = shiftedResverseBitAndEsAndSignificand.range(2+S_WES+N-1-2-S_WES-1, 2+S_WES+N-1-2-S_WES-1-(N-1)+1);
	ap_uint<N> normalOutput = sign.concat(trunctatedEsAndSignificand);
	ap_uint<N-1> zero = 0;
	ap_uint<1> isNaRBit = positValue.getIsNaR();
	ap_uint<N> specialCasesValue = isNaRBit.concat(zero);
	return ((not(positValue.getSignBit()) and not(positValue.getImplicitBit()) or isNaRBit) ? specialCasesValue : normalOutput);
}
