#pragma once
#include "posit_dim.hpp"

#define S_WF PositDim<N>::WF
#define S_WE PositDim<N>::WE
#define S_WES PositDim<N>::WES
#define K_SIZE (S_WE-S_WES)

template<int N>
PositEncoding<N> posit_encoder(PositValue<N> positValue)
{
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

	ap_uint<2> guardBitAndSticky = positValue.getGuardBit().concat(positValue.getStickyBit());
	ap_int<2+S_WES+N-1-2-S_WES+2> reverseBitAndEsAndSignificandAndGuardBitAndSticky = reverseBitAndEsAndSignificand.concat(guardBitAndSticky);
	ap_uint<(1<<(K_SIZE-1))> zeros = 0;
	ap_int<2+S_WES+N-1-2-S_WES+2 + (1<<(K_SIZE-1))> readyToShift = reverseBitAndEsAndSignificandAndGuardBitAndSticky.concat(zeros);

	// printApInt(readyToShift);

	ap_uint<2+S_WES+N-1-2-S_WES+2+ (1<<(K_SIZE-1))> shifted = readyToShift >> absK;

	// printApUint(shifted);

	ap_uint<2+S_WES+N-1-2-S_WES> unroundedResult =  shifted.range(2+S_WES+N-1-2-S_WES+2+ (1<<(K_SIZE-1))-1,2+ (1<<(K_SIZE-1)));

	// printApUint(unroundedResult);

	ap_uint<1> guard = shifted[2+ (1<<(K_SIZE-1))-1];
	ap_uint<1> sticky = not (shifted.range(2+(1<<(K_SIZE-1))-1-1,0) == 0);

	// printApUint(guard);
	// printApUint(sticky);


	ap_uint<1> roundingBit = (guard and not(sticky) and unroundedResult[0]) or (guard and sticky);

	ap_uint<2+S_WES+N-1-2-S_WES> roundedResult = unroundedResult + roundingBit;

	ap_uint<N> normalOutput = sign.concat(roundedResult);
	ap_uint<N-1> zero = 0;
	ap_uint<1> isNaRBit = positValue.getIsNaR();
	ap_uint<N> specialCasesValue = isNaRBit.concat(zero);
	return (((not positValue.getSignBit()) and (not positValue.getImplicitBit())) or isNaRBit) ? specialCasesValue : normalOutput;
}

extern template PositEncoding<8> posit_encoder(PositValue<8> positValue);
extern template PositEncoding<16> posit_encoder(PositValue<16> positValue);
extern template PositEncoding<32> posit_encoder(PositValue<32> positValue);
extern template PositEncoding<64> posit_encoder(PositValue<64> positValue);
