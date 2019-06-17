#include "kulisch_dim.hpp"
#include "tools/static_math.hpp"

using hint::Static_Val;

template<int N>
ap_uint<N> acc_IEEE_rounding(
		KulischAcc<N> acc
){
	#pragma HLS INLINE

	static constexpr int overflow_bits = FPDim<N>::PROD_FP_SPREAD - FPDim<N>::FP_SPREAD;
	static constexpr int underflow_bits = FPDim<N>::ACC_MID - FPDim<N>::FP_SPREAD -FPDim<N>::WF;
	static constexpr int reduced_acc_size = 2*FPDim<N>::FP_SPREAD + FPDim<N>::WF;
	
	static constexpr int ceil_log2_acc_size_minus_1 = Static_Val<reduced_acc_size>::_rlog2-1; 
	static constexpr int pow2_ceil_log2_acc_size_minus_1 = (1<<ceil_log2_acc_size_minus_1); 
	static constexpr int remaining_bits = reduced_acc_size - (1 << (Static_Val<reduced_acc_size>::_rlog2)-1); 
	static constexpr int ceil_log2_remaining_bits = Static_Val<remaining_bits>::_rlog2; 
	static constexpr int pow2_ceil_log2_remaining_bits = (1<<ceil_log2_remaining_bits); 
	
	ap_uint<1> r_s;
	ap_uint<FPDim<N>::WE> r_e;
	ap_uint<FPDim<N>::WF> r_m;


	ap_uint<overflow_bits> high_acc = acc.range(FPDim<N>::ACC_SIZE -1, FPDim<N>::ACC_MID + FPDim<N>::FP_SPREAD);
	ap_uint<reduced_acc_size> mid_acc = acc.range(FPDim<N>::ACC_MID + FPDim<N>::FP_SPREAD -1 , FPDim<N>::ACC_MID - FPDim<N>::FP_SPREAD -FPDim<N>::WF);
	ap_uint<underflow_bits> low_acc = acc.range(FPDim<N>::ACC_MID - FPDim<N>::FP_SPREAD -FPDim<N>::WF -1 , 0);




	r_s = acc[FPDim<N>::ACC_SIZE-1];

	ap_uint<1> overflow_when_neg = high_acc xor r_s;
	ap_uint<1> sticky_low = low_acc.or_reduce();

	ap_uint<Static_Val<reduced_acc_size>::_rlog2 + reduced_acc_size > lzoc_shift = generic_lzoc_shifter<reduced_acc_size>(mid_acc, r_s) ;
	ap_uint<Static_Val<reduced_acc_size>::_rlog2> lzoc = lzoc_shift.range(Static_Val<reduced_acc_size>::_rlog2 + reduced_acc_size-1, reduced_acc_size);
	ap_uint<reduced_acc_size> shifted = lzoc_shift.range(reduced_acc_size-1, 0);
	ap_uint<FPDim<N>::WF> r_m_signed;

	ap_uint<(FPDim<N>::ACC_SIZE - FPDim<N>::WF)> sticky_bits = shifted.range(reduced_acc_size -1-FPDim<N>::WF +1-1,
								reduced_acc_size -1-FPDim<N>::WF +1-1 - (reduced_acc_size - FPDim<N>::WF) +1
		                        );

	ap_uint<1> sticky_tmp = sticky_bits.or_reduce() or sticky_low;

    ap_uint<1> guard1 = shifted[reduced_acc_size -1 -FPDim<N>::WF +1-1-1];
    ap_uint<1> guard2 = shifted[reduced_acc_size -1 -FPDim<N>::WF +1-1];

    ap_uint<1> guard, sticky;

	if(lzoc>(FPDim<N>::SUBNORMAL_LIMIT+1))	{
		r_e = 0;
		r_m_signed = shifted.range(reduced_acc_size -1, reduced_acc_size -1 -FPDim<N>::WF +1) >> (lzoc-(FPDim<N>::SUBNORMAL_LIMIT)-1-1);
		guard=guard2;
		sticky=guard1 or sticky_tmp;

	}
	else{
		r_m_signed = shifted.range(reduced_acc_size -1, reduced_acc_size -1 -FPDim<N>::WF +1-1);
		r_e = FPDim<N>::BIAS-(lzoc)+2 + overflow_bits;
		guard = guard1;
		sticky = sticky_tmp;
	}

	ap_uint<FPDim<N>::WF+1> r_m_rounded;
	if((guard and not(sticky) and not(r_m_signed[0])) or (guard and sticky)){
		r_m_rounded = r_m_signed+1;
	}
	else{
		r_m_rounded = r_m_signed;	
	}

	ap_uint<FPDim<N>::WF+1> r_m_rounded_cut = r_m_rounded.range(FPDim<N>::WF, 0);
	ap_uint<1> overflow = not(r_m_rounded.or_reduce()) and r_s;
	if(r_s){
		r_m = ~r_m_rounded_cut+1;
	}
	else{
		r_m = r_m_rounded_cut;
	}
	r_e += (r_m_rounded[FPDim<N>::WF] and not(r_s)) or overflow;
	ap_uint<1+FPDim<N>::WE> signed_exp = r_s.concat(r_e);


	return signed_exp.concat(r_m);

}


template<int N>
ap_uint<N> signed_acc_IEEE_rounding(
		SignedKulischAcc<N> acc
){
	#pragma HLS INLINE



	static constexpr int overflow_bits = FPDim<N>::PROD_FP_SPREAD - FPDim<N>::FP_SPREAD;
	static constexpr int underflow_bits = FPDim<N>::ACC_MID - FPDim<N>::FP_SPREAD -FPDim<N>::WF;
	static constexpr int reduced_acc_size = 2*FPDim<N>::FP_SPREAD + FPDim<N>::WF;
	
	static constexpr int ceil_log2_acc_size_minus_1 = Static_Val<reduced_acc_size>::_rlog2-1; 
	static constexpr int pow2_ceil_log2_acc_size_minus_1 = (1<<ceil_log2_acc_size_minus_1); 
	static constexpr int remaining_bits = reduced_acc_size - (1 << (Static_Val<reduced_acc_size>::_rlog2)-1); 
	static constexpr int ceil_log2_remaining_bits = Static_Val<remaining_bits>::_rlog2; 
	static constexpr int pow2_ceil_log2_remaining_bits = (1<<ceil_log2_remaining_bits); 
	
	ap_uint<1> r_s;
	ap_uint<FPDim<N>::WE> r_e;
	ap_uint<FPDim<N>::WF> r_m;


	ap_uint<overflow_bits> high_acc = acc.range(FPDim<N>::ACC_SIZE -1, FPDim<N>::ACC_MID + FPDim<N>::FP_SPREAD);
	ap_uint<reduced_acc_size> mid_acc = acc.range(FPDim<N>::ACC_MID + FPDim<N>::FP_SPREAD -1 , FPDim<N>::ACC_MID - FPDim<N>::FP_SPREAD -FPDim<N>::WF);
	ap_uint<underflow_bits> low_acc = acc.range(FPDim<N>::ACC_MID - FPDim<N>::FP_SPREAD -FPDim<N>::WF -1 , 0);




	r_s = acc[FPDim<N>::ACC_SIZE];

	ap_uint<1> overflow_when_neg = high_acc.or_reduce();
	ap_uint<1> sticky_low = low_acc.or_reduce();

	ap_uint<Static_Val<reduced_acc_size>::_rlog2 + reduced_acc_size > lzoc_shift = generic_lzoc_shifter<reduced_acc_size>(mid_acc, 0) ;
	ap_uint<Static_Val<reduced_acc_size>::_rlog2> lzoc = lzoc_shift.range(Static_Val<reduced_acc_size>::_rlog2 + reduced_acc_size-1, reduced_acc_size);
	ap_uint<reduced_acc_size> shifted = lzoc_shift.range(reduced_acc_size-1, 0);
	ap_uint<FPDim<N>::WF> r_m_signed;

	ap_uint<(FPDim<N>::ACC_SIZE - FPDim<N>::WF)> sticky_bits = shifted.range(reduced_acc_size -1-FPDim<N>::WF +1-1,
								reduced_acc_size -1-FPDim<N>::WF +1-1 - (reduced_acc_size - FPDim<N>::WF) +1
		                        );

	ap_uint<1> sticky_tmp = sticky_bits.or_reduce() or sticky_low;

    ap_uint<1> guard1 = shifted[reduced_acc_size -1 -FPDim<N>::WF +1-1-1];
    ap_uint<1> guard2 = shifted[reduced_acc_size -1 -FPDim<N>::WF +1-1];

    ap_uint<1> guard, sticky;

	if(lzoc>(FPDim<N>::SUBNORMAL_LIMIT+1))	{
		r_e = 0;
		r_m_signed = shifted.range(reduced_acc_size -1, reduced_acc_size -1 -FPDim<N>::WF +1) >> (lzoc-(FPDim<N>::SUBNORMAL_LIMIT)-1-1);
		guard=guard2;
		sticky=guard1 or sticky_tmp;

	}
	else{
		r_m_signed = shifted.range(reduced_acc_size -1, reduced_acc_size -1 -FPDim<N>::WF +1-1);
		r_e = FPDim<N>::BIAS-(lzoc)+2 + overflow_bits;
		guard = guard1;
		sticky = sticky_tmp;
	}

	ap_uint<FPDim<N>::WF+1> r_m_rounded;
	if((guard and not(sticky) and not(r_m_signed[0])) or (guard and sticky)){
		r_m_rounded = r_m_signed+1;
	}
	else{
		r_m_rounded = r_m_signed;	
	}

	ap_uint<FPDim<N>::WF+1> r_m_rounded_cut = r_m_rounded.range(FPDim<N>::WF, 0);
	// ap_uint<1> overflow = not(r_m_rounded.or_reduce()) and r_s;
	// if(r_s){
	// 	r_m = ~r_m_rounded_cut+1;
	// }
	// else{
		r_m = r_m_rounded_cut;
	// }
	r_e += r_m_rounded[FPDim<N>::WF];
	ap_uint<1+FPDim<N>::WE> signed_exp = r_s.concat(r_e);


	return signed_exp.concat(r_m);


}
