#pragma once
#include <cstdio>

#include "ap_int.h"
#include "shifter.hpp"
#include "posit_dim.hpp"


template<int N>
Quire<N> add_sub_quire(
		Quire<N> quire, 
		PositProd<N> input,
		ap_uint<1> isSub
){
	#pragma HLS INLINE
	static constexpr int LOG2_EXT_SUM_SIZE = Static_Val<PositDim<N>::ExtQuireSize-1 +PositDim<N>::ProdSignificandSize>::_log2;

	ap_int<PositDim<N>::ProdSignificandSize+1> inputSignificand = input.getSignedSignificand();
	ap_int<PositDim<N>::ProdSignificandSize+1> complementedInputIfIsSub;
	ap_uint<1> sign = inputSignificand[PositDim<N>::ProdSignificandSize];

	for (int i=0; i<PositDim<N>::ProdSignificandSize+1; i++){
	#pragma HLS UNROLL
		complementedInputIfIsSub[i] = input[i] ^ isSub;
	}

	ap_uint<PositDim<N>::ProdExpSize> shiftValue = input.getExp();
	ap_int<1<<LOG2_EXT_SUM_SIZE> ext = complementedInputIfIsSub;

	ap_int<(1<<LOG2_EXT_SUM_SIZE)> shiftedInput = shifter<LOG2_EXT_SUM_SIZE>(ext,shiftValue,isSub);
	ap_int<PositDim<N>::ExtQuireSize-1> shiftedInputShrinked = shiftedInput.range(PositDim<N>::ExtQuireSize-1 + PositDim<N>::ProdSignificandSize-1, PositDim<N>::ProdSignificandSize);

	ap_uint<PositDim<N>::ExtQuireSize-1> quireWithoutSignAndNARBit = quire.getQuireWithoutNaR();

	ap_uint<PositDim<N>::ExtQuireSize-1> sumResult = shiftedInputShrinked + quireWithoutSignAndNARBit + isSub;

	ap_uint<1> resultIsNaR = quire.getIsNaR() || input.getIsNaR();

	return Quire<N>(resultIsNaR.concat(sumResult));
}

extern template Quire<8> add_sub_quire<8>(Quire<8>, PositProd<8>, ap_uint<1>);
extern template Quire<16> add_sub_quire<16>(Quire<16>, PositProd<16>, ap_uint<1>);
extern template Quire<32> add_sub_quire<32>(Quire<32>, PositProd<32>, ap_uint<1>);
extern template Quire<64> add_sub_quire<64>(Quire<64>, PositProd<64>, ap_uint<1>);

template<int bankSize> 
constexpr int getIndex(int index, bool isUpper)
{
	return bankSize*index - isUpper;
}

template <int N, int bankSize, int spread>
ap_uint<bankSize> getToAddRec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<PositDim<N>::ProdExpSize - getShiftSize<bankSize>()> stageSelect,
    ap_uint<1> inputSign,
    ap_uint<1> isSub,
    ap_int<getExtShiftSize<N, bankSize>()> shiftedSignificand,
    typename enable_if<(spread < 1)>::type* dummy = 0
)
{
	#pragma HLS INLINE
	
	static constexpr int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositDim<N>::ProdSignificandSize,bankSize>::val;

	if((stageIndex>=(stageSelect+getMantSpread<N, bankSize>()-BANKS_FOR_USELESS_BITS) && (inputSign ^ isSub)) || ((stageIndex<stageSelect) && isSub )){
		return -1;
	}
	else{
		return 0;
	}
}

template <int N, int bankSize, int spread>
ap_uint<bankSize> getToAddRec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<PositDim<N>::ProdExpSize - getShiftSize<bankSize>()> stageSelect,
    ap_uint<1> inputSign,
    ap_uint<1> isSub,
    ap_int<getExtShiftSize<N, bankSize>()> shiftedSignificand,
    typename enable_if<(spread >= 1)>::type* dummy = 0
    )
{	
	#pragma HLS INLINE
	static constexpr int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositDim<N>::ProdSignificandSize,bankSize>::val;

    if (stageIndex == (stageSelect+spread-1-BANKS_FOR_USELESS_BITS)) {
        return shiftedSignificand.range((spread*bankSize)-1,(spread-1)*bankSize);
    } else {
        return getToAddRec<N, bankSize, (spread-1)>(stageIndex, stageSelect, inputSign, isSub, shiftedSignificand);
    }
}


template <int N, int bankSize, int spread>
ap_uint<bankSize> getToAdd(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<PositDim<N>::ProdExpSize - getShiftSize<bankSize>()> stageSelect,
    ap_uint<1> inputSign,
    ap_uint<1> isSub,
    ap_int<getExtShiftSize<N, bankSize>()> shiftedSignificand
		){
	#pragma HLS INLINE
	return getToAddRec<N, bankSize, spread>(stageIndex, stageSelect, inputSign, isSub, shiftedSignificand);
}

template<int N, int bankSize>
ap_uint<bankSize+1> add_sub_quire_stage(SegmentedQuire<N, bankSize> quire, 
										ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
										ap_uint<PositDim<N>::ProdExpSize - getShiftSize<bankSize>()> stageSelect,
										ap_uint<1> inputSign, 
										ap_uint<1> isSub,
										ap_int<getExtShiftSize<N, bankSize>()> shiftedSignificand)
{
	#pragma HLS INLINE
	ap_uint<bankSize+1> quireBank = quire.getBank(stageIndex);
	
	ap_uint<1> quireCarry = (stageIndex==0) ? isSub : quire.getCarry(stageIndex-1);

	ap_uint<bankSize> toAdd = getToAdd<N, bankSize, getMantSpread<N, bankSize>()>(stageIndex, stageSelect, inputSign, isSub, shiftedSignificand);

	ap_uint<bankSize+1> sum = quireBank + toAdd + quireCarry;
	return sum;
}


template<int N, int bankSize>
SegmentedQuire<N, bankSize> segmented_add_sub_quire(SegmentedQuire<N, bankSize> quire, 
													PositProd<N> input,
													ap_uint<1> isSub)
{	
	#pragma HLS INLINE

	static constexpr int SHIFTED_SIGNIFICAND_SIZE = PositDim<N>::ProdSignificandSize+1 + (1<<getShiftSize<bankSize>());
	static constexpr int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositDim<N>::ProdSignificandSize,bankSize>::val;
	static constexpr int ENCODING_BITS_FOR_USELESS_BANKS = Static_Val<BANKS_FOR_USELESS_BITS>::_log2;
	static constexpr int padding = bankSize*BANKS_FOR_USELESS_BITS - PositDim<N>::ProdSignificandSize;
	static constexpr int LOG2_SHIFT_SIZE = Static_Val<getExtShiftSize<N, bankSize>()>::_log2;

	ap_uint<PositDim<N>::ProdExpSize+ENCODING_BITS_FOR_USELESS_BANKS> prodExp = input.getExp()+padding;
	ap_int<PositDim<N>::ProdSignificandSize+1> inputSignificand = input.getSignedSignificand();
	ap_uint<1> sign = inputSignificand[PositDim<N>::ProdSignificandSize];
	
	ap_int<PositDim<N>::ProdSignificandSize+1> complementedInputIfIsSub;
	for (int i=0; i<PositDim<N>::ProdSignificandSize+1; i++){
		#pragma HLS UNROLL
		complementedInputIfIsSub[i] = input[i] ^ isSub;
	}

	ap_uint<getShiftSize<bankSize>()> shiftValue = prodExp.range(getShiftSize<bankSize>()-1,0);

	ap_int<1<<LOG2_SHIFT_SIZE> ext = complementedInputIfIsSub;
	ap_int<(1<<LOG2_SHIFT_SIZE)> shiftedInput = shifter<LOG2_SHIFT_SIZE>(ext,shiftValue,isSub);
	ap_int<getExtShiftSize<N, bankSize>()> shiftedInputShrinked = shiftedInput.range(getExtShiftSize<N, bankSize>()-1, 0);

	ap_uint<PositDim<N>::ProdExpSize - getShiftSize<bankSize>()> stageSelect = prodExp.range(PositDim<N>::ProdExpSize-1, getShiftSize<bankSize>());

	SegmentedQuire<N, bankSize> fullQuire = SegmentedQuire<N, bankSize>(0);

	for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
		#pragma HLS UNROLL
		ap_uint<bankSize+1> stageResult = add_sub_quire_stage<N,bankSize>(quire, i, stageSelect, sign, isSub, shiftedInputShrinked);
		fullQuire[i] = stageResult[bankSize];
		fullQuire.range(getIndex<bankSize>(i+1, 1)+getNbStages<N, bankSize>(), getIndex<bankSize>(i, 0)+getNbStages<N, bankSize>()) = stageResult.range(bankSize-1,0);
	}
	ap_uint<1> resultIsNaR = input.getIsNaR() || quire.getIsNaR();
	ap_uint<PositDim<N>::ExtQuireSize+getNbStages<N, bankSize>()-1> fullQuireWoNaR = fullQuire.range(PositDim<N>::ExtQuireSize+getNbStages<N, bankSize>()-1-1, 0);
	return SegmentedQuire<N, bankSize>(resultIsNaR.concat(fullQuireWoNaR));
}

template<int N, int bankSize>
Quire<N> propagateCarries(SegmentedQuire<N, bankSize> quire)
{	
	#pragma HLS INLINE
	SegmentedQuire<N, bankSize> fullQuire = quire;
	for(int j=0; j<getNbStages<N, bankSize>(); j++){
		for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
			#pragma HLS UNROLL
			ap_uint<bankSize+1> stageResult = add_sub_quire_stage<N,bankSize>(fullQuire, i, 0, 0, 0, 0);
			fullQuire[i] = stageResult[bankSize];
			fullQuire.range(getIndex<bankSize>(i+1, 1)+getNbStages<N, bankSize>(), getIndex<bankSize>(i, 0)+getNbStages<N, bankSize>()) = stageResult.range(bankSize-1,0);
		}	
	}
	return fullQuire.getAsQuireWoCarries();
}
