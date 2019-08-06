#pragma once
#include <cstdio>

#include "ap_int.h"
#include "posit_dim.hpp"
#include "bitvector/lzoc_shifter.hpp"

template<int N, int WES, int NB_CARRY>
Quire<N, WES, NB_CARRY> add_sub_quire(
        Quire<N, WES, NB_CARRY> quire,
        PositProd<N, WES> input,
		ap_uint<1> isSub
){
	#pragma HLS INLINE
    static constexpr int LOG2_EXT_SUM_SIZE = Static_Val<quire.Size-1 +input.SignificandSize>::_log2;

    ap_int<input.SignificandSize+1> inputSignificand = input.getSignedSignificand();
    ap_int<input.SignificandSize+1> complementedInputIfIsSub;
    ap_uint<1> sign = inputSignificand[input.SignificandSize];

    for (int i=0; i<input.SignificandSize + 1; i++){
	#pragma HLS UNROLL
		complementedInputIfIsSub[i] = input[i] ^ isSub;
	}

    ap_uint<input.ExpSize> shiftValue = input.getExp();
	ap_int<1<<LOG2_EXT_SUM_SIZE> ext = complementedInputIfIsSub;

	ap_int<(1<<LOG2_EXT_SUM_SIZE)> shiftedInput = shifter<LOG2_EXT_SUM_SIZE>(ext,shiftValue,isSub);
    ap_int<quire.Size-1> shiftedInputShrinked = shiftedInput.range(quire.Size + input.SignificandSize-2, input.SignificandSize);

    ap_uint<quire.Size - 1> quireWithoutSignAndNARBit = quire.getQuireWithoutNaR();

    ap_uint<quire.Size-1> sumResult = shiftedInputShrinked + quireWithoutSignAndNARBit + isSub;

    ap_uint<1> resultIsNaR = quire.getIsNaR() | input.getIsNaR();

    return Quire<N, WES, NB_CARRY>{resultIsNaR.concat(sumResult)};
}


template<int bankSize> 
constexpr int getIndex(int index, bool isUpper)
{
	return bankSize*index - isUpper;
}

template <int N, int WES, int NB_CARRY, int bankSize, int spread>
ap_uint<bankSize> getToAddRec(
    ap_uint<Static_Val<getNbStages<N, WES, NB_CARRY, bankSize>()>::_log2> stageIndex,
    ap_uint<PositProd<N, WES>::ExpSize - getShiftSize<bankSize>()> stageSelect,
    ap_uint<1> inputSign,
    ap_uint<1> isSub,
    ap_int<getExtShiftSize<N, WES, bankSize>()>,
    typename enable_if<(spread < 1)>::type* = 0
)
{
	#pragma HLS INLINE
	
    static constexpr int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositProd<N, WES>::SignificandSize, bankSize>::val;

    if((stageIndex>=(stageSelect+getMantSpread<N, WES, bankSize>()-BANKS_FOR_USELESS_BITS) && (inputSign ^ isSub)) || ((stageIndex<stageSelect) && isSub )){
		return -1;
	}
	else{
		return 0;
	}
}

template <int N, int WES, int NB_CARRY, int bankSize, int spread>
ap_uint<bankSize> getToAddRec(
    ap_uint<Static_Val<getNbStages<N, WES, NB_CARRY, bankSize>()>::_log2> stageIndex,
    ap_uint<PositProd<N, WES>::ExpSize - getShiftSize<bankSize>()> stageSelect,
    ap_uint<1> inputSign,
    ap_uint<1> isSub,
    ap_int<getExtShiftSize<N, WES, bankSize>()> shiftedSignificand,
    typename enable_if<(spread >= 1)>::type* = 0
    )
{	
	#pragma HLS INLINE
    static constexpr int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositProd<N, WES>::SignificandSize, bankSize>::val;

    if (stageIndex == (stageSelect+spread-1-BANKS_FOR_USELESS_BITS)) {
        return shiftedSignificand.range((spread*bankSize)-1,(spread-1)*bankSize);
    } else {
        return getToAddRec<N, WES, NB_CARRY, bankSize, (spread-1)>(stageIndex, stageSelect, inputSign, isSub, shiftedSignificand);
    }
}


template <int N, int WES, int NB_CARRY, int bankSize, int spread>
ap_uint<bankSize> getToAdd(
    ap_uint<Static_Val<getNbStages<N, WES, NB_CARRY, bankSize>()>::_log2> stageIndex,
    ap_uint<PositProd<N, WES>::ExpSize - getShiftSize<bankSize>()> stageSelect,
    ap_uint<1> inputSign,
    ap_uint<1> isSub,
    ap_int<getExtShiftSize<N, WES, bankSize>()> shiftedSignificand
		){
	#pragma HLS INLINE
    return getToAddRec<N, WES, NB_CARRY, bankSize, spread>(stageIndex, stageSelect, inputSign, isSub, shiftedSignificand);
}

template<int N, int WES, int NB_CARRY, int bankSize>
ap_uint<bankSize+1> add_sub_quire_stage(SegmentedQuire<N, WES, NB_CARRY, bankSize> quire,
                                        ap_uint<Static_Val<getNbStages<N, WES, NB_CARRY, bankSize>()>::_log2> stageIndex,
                                        ap_uint<PositProd<N, WES>::ExpSize - getShiftSize<bankSize>()> stageSelect,
										ap_uint<1> inputSign, 
										ap_uint<1> isSub,
                                        ap_int<getExtShiftSize<N, WES, bankSize>()> shiftedSignificand)
{
	#pragma HLS INLINE
	ap_uint<bankSize+1> quireBank = quire.getBank(stageIndex);
	
	ap_uint<1> quireCarry = (stageIndex==0) ? isSub : quire.getCarry(stageIndex-1);

    ap_uint<bankSize> toAdd = getToAdd<N, WES, NB_CARRY, bankSize, getMantSpread<N, WES, bankSize>()>(stageIndex, stageSelect, inputSign, isSub, shiftedSignificand);

	ap_uint<bankSize+1> sum = quireBank + toAdd + quireCarry;
	return sum;
}

template<int NB_CARRY, int bankSize, int N, int WES>
SegmentedQuire<N, WES, NB_CARRY, bankSize> segmented_add_sub_quire(SegmentedQuire<N, WES, NB_CARRY, bankSize> quire,
                                                    PositProd<N, WES> input,
													ap_uint<1> isSub)
{	
	#pragma HLS INLINE

    //static constexpr int SHIFTED_SIGNIFICAND_SIZE = PositProd<N, WES>::SignificandSize+1 + (1<<getShiftSize<bankSize>());
    static constexpr int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositProd<N, WES>::SignificandSize, bankSize>::val;
	static constexpr int ENCODING_BITS_FOR_USELESS_BANKS = Static_Val<BANKS_FOR_USELESS_BITS>::_log2;
    static constexpr int padding = bankSize*BANKS_FOR_USELESS_BITS - PositProd<N, WES>::SignificandSize;
    static constexpr int LOG2_SHIFT_SIZE = Static_Val<getExtShiftSize<N, WES, bankSize>()>::_log2;
    static constexpr int NB_STAGES = getNbStages<N, WES, NB_CARRY, bankSize>();

    ap_uint<PositProd<N, WES>::ExpSize+ENCODING_BITS_FOR_USELESS_BANKS> prodExp = input.getExp()+padding;
    ap_int<PositProd<N, WES>::SignificandSize+1> inputSignificand = input.getSignedSignificand();
    ap_uint<1> sign = inputSignificand[PositProd<N, WES>::SignificandSize];
	
    ap_int<PositProd<N, WES>::SignificandSize+1> complementedInputIfIsSub;
    for (int i=0; i<PositProd<N, WES>::SignificandSize + 1; i++){
		#pragma HLS UNROLL
		complementedInputIfIsSub[i] = input[i] ^ isSub;
	}

	ap_uint<getShiftSize<bankSize>()> shiftValue = prodExp.range(getShiftSize<bankSize>()-1,0);

	ap_int<1<<LOG2_SHIFT_SIZE> ext = complementedInputIfIsSub;
	ap_int<(1<<LOG2_SHIFT_SIZE)> shiftedInput = shifter<LOG2_SHIFT_SIZE>(ext,shiftValue,isSub);
    ap_int<getExtShiftSize<N, WES, bankSize>()> shiftedInputShrinked = shiftedInput.range(getExtShiftSize<N, WES, bankSize>()-1, 0);

    ap_uint<PositProd<N, WES>::ExpSize - getShiftSize<bankSize>()> stageSelect = prodExp.range(PositProd<N, WES>::ExpSize-1, getShiftSize<bankSize>());

    SegmentedQuire<N, WES, NB_CARRY, bankSize> fullQuire{0};

    for(int i=NB_STAGES-1; i>=0; i--){
		#pragma HLS UNROLL
        ap_uint<bankSize+1> stageResult = add_sub_quire_stage<N, WES, NB_CARRY, bankSize>(
                    quire,
                    i,
                    stageSelect,
                    sign,
                    isSub,
                    shiftedInputShrinked
            );
		fullQuire[i] = stageResult[bankSize];
        fullQuire.range(getIndex<bankSize>(i+1, 1)+NB_STAGES, getIndex<bankSize>(i, 0)+NB_STAGES) = stageResult.range(bankSize-1,0);
	}
	ap_uint<1> resultIsNaR = input.getIsNaR() || quire.getIsNaR();
    ap_uint<SegmentedQuire<N, WES, NB_CARRY, bankSize>::Size - 1> fullQuireWoNaR = fullQuire.range(Quire<N, WES, NB_CARRY>::Size+NB_STAGES-2, 0);
    return SegmentedQuire<N, WES, NB_CARRY, bankSize>(resultIsNaR.concat(fullQuireWoNaR));
}

template<int N, int WES, int NB_CARRY, int bankSize>
Quire<N, WES, NB_CARRY> propagateCarries(SegmentedQuire<N, WES, NB_CARRY, bankSize> quire)
{	
	#pragma HLS INLINE
    constexpr int NB_STAGES = getNbStages<N, WES, NB_CARRY, bankSize>();
    SegmentedQuire<N, WES, NB_CARRY, bankSize> fullQuire = quire;
    for(int j=0; j<NB_STAGES; j++){
        for(int i=NB_STAGES-1; i>=0; i--){
			#pragma HLS UNROLL
            ap_uint<bankSize+1> stageResult = add_sub_quire_stage<N, WES, NB_CARRY, bankSize>(
                        fullQuire,
                        ap_uint<Static_Val<NB_STAGES>::_log2>{i},
                        ap_uint<PositProd<N, WES>::ExpSize - getShiftSize<bankSize>()>{0},
                        ap_uint<1>{0},
                        ap_uint<1>{0},
                        ap_int<getExtShiftSize<N, WES, bankSize>()>{0}
                );
			fullQuire[i] = stageResult[bankSize];
            fullQuire.range(getIndex<bankSize>(i+1, 1)+NB_STAGES,
                            getIndex<bankSize>(i, 0)+NB_STAGES) = stageResult.range(bankSize-1,0);
		}	
	}
	return fullQuire.getAsQuireWoCarries();
}
