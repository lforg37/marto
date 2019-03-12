#include "kulisch_dim.hpp"
#include "lzoc_shifter.hpp"
#include "utils.hpp"
#include "shifter.hpp"
#include <cstdio>
using namespace std;

template<int N>
FPProd<N> exact_prod(ap_uint<N> in1,
		ap_uint<N> in2
){
	#pragma HLS INLINE
	ap_uint<FPDim<N>::WF+1> m1, m2;
	ap_uint<FPDim<N>::WE> e1, e2;
	ap_uint<1> s1, s2;
	ap_uint<2*FPDim<N>::WF+2> mult_m;
	ap_uint<FPDim<N>::WE+1> mult_e;

	m1 = in1.range(FPDim<N>::WF-1, 0);
	e1 = in1.range(FPDim<N>::WF+FPDim<N>::WE-1, FPDim<N>::WF);
	m2 = in2.range(FPDim<N>::WF-1, 0);
	e2 = in2.range(FPDim<N>::WF+FPDim<N>::WE-1, FPDim<N>::WF);
	s1 = in1[FPDim<N>::WF+FPDim<N>::WE];
	s2 = in2[FPDim<N>::WF+FPDim<N>::WE];

	if(e1 == 0){
		m1[FPDim<N>::WF] = 0;
	}
	else{
		m1[FPDim<N>::WF] = 1;
	}

	if(e2 == 0){
		m2[FPDim<N>::WF] = 0;
	}
	else{
		m2[FPDim<N>::WF] = 1;
	}

	ap_uint<1> mult_s = s1 ^ s2;

	mult_e = (ap_uint<FPDim<N>::WE+1>)e1 + (ap_uint<FPDim<N>::WE+1>)e2;
	mult_m = ((ap_uint<2*FPDim<N>::WF+2>) m1) * ((ap_uint<2*FPDim<N>::WF+2>) m2);
	return FPProd<N>(	mult_s,
						mult_e,
						mult_m	
					);
}


template<int N>
KulischAcc<N> kulisch_accumulator(
		KulischAcc<N> acc, 
		FPProd<N> prod
){
	#pragma HLS INLINE

	// cout << "Prod :\n";
	// prod.printContent();

	ap_int<2*FPDim<N>::WF+3> ext_significand = prod.getSignificand();

	if(prod.getSignBit() == 1){
		ext_significand = ~((ap_int<2*FPDim<N>::WF+3>)prod.getSignificand()) + 1;
	}
	else{
		ext_significand = prod.getSignificand();
	}

	ap_uint<FPDim<N>::ACC_SIZE> shifted = ((ap_uint<FPDim<N>::ACC_SIZE>) ext_significand) << (prod.getExp());
	ap_uint<FPDim<N>::ACC_SIZE> r_acc = shifted + acc;

	// cout << "Acc " << FPDim<N>::ACC_SIZE << ":\n";
	// printApUint(r_acc);

	return r_acc;
}

template<int N>
ap_uint<N> acc_to_fp(
		KulischAcc<N> acc
){
	#pragma HLS INLINE

	static constexpr int ceil_log2_acc_size = Static_Val<FPDim<N>::ACC_SIZE>::_log2; 
	static constexpr int added_bits = (1<<ceil_log2_acc_size) - FPDim<N>::ACC_SIZE;
	ap_uint<1> r_s;
	ap_uint<FPDim<N>::WE> r_e;
	ap_uint<FPDim<N>::WF> r_m;


	ap_int<(1<<ceil_log2_acc_size)> ext_acc = (ap_int<FPDim<N>::ACC_SIZE>)acc;
	// fprintf(stderr, "ext acc %d\n", (int)ext_acc);
	// printApInt(ext_acc);
	r_s = acc[FPDim<N>::ACC_SIZE-1];

	ap_uint<ceil_log2_acc_size + (1<<ceil_log2_acc_size)> lzoc_shift = lzoc_shifter<ceil_log2_acc_size>(ext_acc, r_s);

	ap_uint<ceil_log2_acc_size> lzoc = lzoc_shift.range(ceil_log2_acc_size + (1<<ceil_log2_acc_size)-1, (1<<ceil_log2_acc_size));
	ap_uint<(1<<ceil_log2_acc_size)> shifted = lzoc_shift.range((1<<ceil_log2_acc_size)-1, 0);
	
	// fprintf(stderr, "lzoc %d\n", (int)lzoc);
	// printApUint(lzoc);
	// fprintf(stderr, "shifted \n");
	// printApUint(shifted);

	// fprintf(stderr, "Limit : %d\n", FPDim<N>::SUBNORMAL_LIMIT+added_bits);
	// fprintf(stderr, "SUBNORMAL_LIMIT : %d\n", FPDim<N>::SUBNORMAL_LIMIT);
	// fprintf(stderr, "added_bits : %d\n", added_bits);
	ap_uint<FPDim<N>::WF> r_m_signed;

	ap_uint<(FPDim<N>::ACC_SIZE - FPDim<N>::WF)> sticky_bits = shifted.range((1<<ceil_log2_acc_size)-1-FPDim<N>::WF +1-1-1,
								(1<<ceil_log2_acc_size)-1-FPDim<N>::WF +1-1-1 - (FPDim<N>::ACC_SIZE - FPDim<N>::WF) +1
		                        );

	ap_uint<1> sticky_tmp = not(sticky_bits== 0);

    ap_uint<1> guard1 = shifted[(1<<ceil_log2_acc_size)-1 -FPDim<N>::WF +1-1-1];
    ap_uint<1> guard2 = shifted[(1<<ceil_log2_acc_size)-1 -FPDim<N>::WF +1-1];

    ap_uint<1> guard, sticky;

	if(lzoc>(FPDim<N>::SUBNORMAL_LIMIT+added_bits+1))	{
		r_e = 0;
		r_m_signed = shifted.range((1<<ceil_log2_acc_size)-1, (1<<ceil_log2_acc_size)-1 -FPDim<N>::WF +1) >> (lzoc-(FPDim<N>::SUBNORMAL_LIMIT+added_bits)-1-1);
		guard=guard2;
		sticky=guard1 or sticky_tmp;
	}
	else{
		r_m_signed = shifted.range((1<<ceil_log2_acc_size)-1-1, (1<<ceil_log2_acc_size)-1 -FPDim<N>::WF +1-1);
		r_e = FPDim<N>::BIAS-(lzoc-added_bits)+2;
		guard = guard1;
		sticky = sticky_tmp;
	}
	// printApUint(r_m_signed);
	// printApUint(sticky_bits);
	// printApUint(shifted);

	// fprintf(stderr, "guard %d, sticky %d\n",  (int)guard, (int)sticky);
	// fprintf(stderr, "guard1 %d, guard2 %d, sticky_tmp %d\n",  (int)guard1, (int)guard2, (int)sticky_tmp);
	ap_uint<FPDim<N>::WF+1> r_m_rounded;
	if((guard and not(sticky) and not(r_m_signed[0])) or (guard and sticky)){
		// fprintf(stderr, "rounded\n");
		r_m_rounded = r_m_signed+1;
			// fprintf(stderr, "round\n");
			// if(r_m_rounded[FPDim<N>::WF]==1){
			// 	fprintf(stderr, "round overflowed\n");
			// }
	}
	else{
		// fprintf(stderr, "not round\n");
		r_m_rounded = r_m_signed;
		// fprintf(stderr, "not rounded\n");		
	}

	ap_uint<FPDim<N>::WF+1> r_m_rounded_cut = r_m_rounded.range(FPDim<N>::WF, 0);
	ap_uint<1> overflow = (r_m_rounded==0) and r_s;
	if(r_s){
		r_m = ~r_m_rounded_cut+1;
	}
	else{
		r_m = r_m_rounded_cut;
	}
	// fprintf(stderr, "r_s %d\n", (int)r_s);
	r_e+=(r_m_rounded[FPDim<N>::WF] and not(r_s)) or overflow;
	ap_uint<1+FPDim<N>::WE> signed_exp = r_s.concat(r_e);


	return signed_exp.concat(r_m);

}












template<int N, int bankSize>
static constexpr int ext_shift_size(){
	return Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val * bankSize;
}

template<int bankSize> 
constexpr int getIndex(int index, bool isUpper)
{
	return bankSize*index - isUpper;
}

template <int N, int bankSize, int spread>
ap_uint<bankSize> getToAddRec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)> shiftedSignificand,
    typename enable_if<(spread < 1)>::type* dummy = 0
)
{
	#pragma HLS INLINE
	
	if((stageIndex>=(stageSelect+getMantSpread<N, bankSize>()) && inputSign)){
		return -1;
	}
	else{
		return 0;
	}
}

template <int N, int bankSize, int spread>
ap_uint<bankSize> getToAddRec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)> shiftedSignificand,
    typename enable_if<(spread >= 1)>::type* dummy = 0
    )
{	
	#pragma HLS INLINE
	// fprintf(stderr, "if (stageIndex (%d) == (stageSelect (%d) + spread (%d)-1) (%d))\n", (int)stageIndex, (int)stageSelect, (int)spread, (int)(stageSelect+spread-1));
    if (stageIndex == (stageSelect+spread-1)) {
    	// fprintf(stderr, "size : %d\n", (int)(1<<Static_Val<spread*bankSize>::_log2));
    	// fprintf(stderr, "input :");
    	// printApInt(shiftedSignificand);
        ap_uint<bankSize> tmp =shiftedSignificand.range(((spread)*bankSize)-1,(spread-1)*bankSize); 
    	// fprintf(stderr, "taking from %d to %d:", (int) ((spread)*bankSize)-1, (int)(spread-1)*bankSize);
    	// printApUint(tmp);
        return tmp;
    } else {
        return getToAddRec<N, bankSize, (spread-1)>(stageIndex, stageSelect, inputSign, shiftedSignificand);
    }
}


template <int N, int bankSize, int spread>
ap_uint<bankSize> getToAdd(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)	> shiftedSignificand
		){
	#pragma HLS INLINE
	return getToAddRec<N, bankSize, spread>(stageIndex, stageSelect, inputSign, shiftedSignificand);
}

template<int N, int bankSize>
ap_uint<bankSize+1> add_sub_acc_stage(SegmentedKulischAcc<N, bankSize> acc, 
										ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
										ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
										ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)> shiftedSignificand,
										ap_uint<1> sign,
										ap_uint<1> carry0 =0)
{
	#pragma HLS INLINE
	ap_uint<bankSize+1> bank = acc.getBank(stageIndex);
	// fprintf(stderr, "main input with size %d: ", (int) (1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2));
	// printApInt(shiftedSignificand);
	
	ap_uint<1> accCarry = (stageIndex==0) ? carry0 : acc.getCarry(stageIndex-1);

	ap_uint<bankSize> toAdd = getToAdd<N, bankSize, getMantSpread<N, bankSize>()>(stageIndex, stageSelect, sign, shiftedSignificand);

	ap_uint<bankSize+1> sum = bank + toAdd + accCarry;
	return sum;
}




template<int N, int bankSize>
SegmentedKulischAcc<N, bankSize> segmented_kulisch_accumulator(
													SegmentedKulischAcc<N, bankSize> acc, 
													FPProd<N> prod
													)
{	
	#pragma HLS INLINE

	// fprintf(stderr, "acc size %d, ext acc size %d, carry %d\n",FPDim<N>::ACC_SIZE, getSegmentedAccSize<N, bankSize>(), getNbStages<N, bankSize>());

	static constexpr int spread = Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val+1;
	static constexpr int bits_for_shift = Static_Val<spread*bankSize>::_log2;
	// fprintf(stderr, "spread : %d\n", spread);

	// fprintf(stderr, "input exp : ");
	// printApUint(prod.getExp());
	// fprintf(stderr, "bias added : %d\n", 2*FPDim<N>::WF+2);
	// ap_uint<FPDim<N>::WE+1 +1> prodExp = (ap_uint<FPDim<N>::WE+1 +1>)(prod.getExp())+2*FPDim<N>::WF+2;
	ap_uint<FPDim<N>::WE+1 +1> prodExp = (ap_uint<FPDim<N>::WE+1 +1>)(prod.getExp());

	// fprintf(stderr, "prod exp : ");
	// printApUint(prodExp);
	// ap_uint<FPDim<N>::WE+1 +1> prodExp = (ap_uint<FPDim<N>::WE+1 +1>)(prod.getExp());
	ap_uint<2*FPDim<N>::WF+2> inputSignificand = prod.getSignificand();
	ap_uint<2*FPDim<N>::WF+3> inputSignificandComplemented;
	if(prod.getSignBit() == 1){
		inputSignificandComplemented = ~((ap_uint<2*FPDim<N>::WF+3>)inputSignificand) + 1;
	}
	else{
		inputSignificandComplemented = inputSignificand;
	}

	ap_int<2*FPDim<N>::WF+2+2> inputSignificandWithSign = prod.getSignBit().concat(inputSignificandComplemented);
	// fprintf(stderr, "inputSignificand : ");
	// printApUint(inputSignificand);
	// fprintf(stderr, "inputSignificandWithSign : ");
	// printApInt(inputSignificandWithSign);

	Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val *bankSize ;


	ap_uint<Static_Val<bankSize>::_log2> shiftValue = prodExp.range(Static_Val<bankSize>::_log2-1,0);
	// fprintf(stderr, "shiftValue : ");
	// printApUint(shiftValue);

	ap_int<(1<<bits_for_shift)> ext = inputSignificandWithSign;
	// fprintf(stderr, "ext size : %d\n", (1<<bits_for_shift));
	// fprintf(stderr, "ext : ");
	// printApInt(ext);

	ap_int<(1<<bits_for_shift)> shiftedInput = shifter<bits_for_shift>(ext,shiftValue,0);

	// fprintf(stderr, "shiftedInput : ");
	// printApInt(shiftedInput);


	ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect = prodExp.range(FPDim<N>::WE+1 +1 -1, Static_Val<bankSize>::_log2);

	// fprintf(stderr, "stageSelect : ");
	// printApUint(stageSelect);
	SegmentedKulischAcc<N, bankSize> fullAcc = SegmentedKulischAcc<N, bankSize>(0, 0);

	for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
		#pragma HLS UNROLL
		ap_uint<bankSize+1> stageResult = add_sub_acc_stage<N,bankSize>(acc, i, stageSelect, shiftedInput, prod.getSignBit());
		fullAcc[i] = stageResult[bankSize];
		// fprintf(stderr, "high : %d, low : %d\n", getIndex<bankSize>(i+1, 1)+getNbStages<N, bankSize>(),  getIndex<bankSize>(i, 0)+getNbStages<N, bankSize>());
		fullAcc.range(getIndex<bankSize>(i+1, 1)+getNbStages<N, bankSize>(), getIndex<bankSize>(i, 0)+getNbStages<N, bankSize>()) = stageResult.range(bankSize-1,0);
	}

	return fullAcc;
}

template<int N, int bankSize>
KulischAcc<N> Kulisch_propagate_carries(SegmentedKulischAcc<N, bankSize> acc)
{	
	#pragma HLS INLINE
	SegmentedKulischAcc<N, bankSize> fullacc = acc;
	for(int j=0; j<getNbStages<N, bankSize>()+1; j++){
		for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
			#pragma HLS UNROLL
			ap_uint<bankSize+1> stageResult = add_sub_acc_stage<N, bankSize>(fullacc, (ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> ) i, (ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1>) 0, (ap_uint<Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val * bankSize>) 0, 0);
			// ap_uint<bankSize+1> stageResult = add_sub_acc_stage<N, bankSize>(fullacc, (ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> ) i, (ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1>) 0, (ap_uint<Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val * bankSize>) 0, (ap_uint<1>)((j==0) and (i==0) and acc.isNeg()));
			fullacc[i] = stageResult[bankSize];
			fullacc.range(getIndex<bankSize>(i+1, 1)+getNbStages<N, bankSize>(), getIndex<bankSize>(i, 0)+getNbStages<N, bankSize>()) = stageResult.range(bankSize-1,0);
		}	
	}
	// fullacc.printContent();
	return fullacc.range(FPDim<N>::ACC_SIZE + getNbStages<N, bankSize>()-1, getNbStages<N, bankSize>());
}