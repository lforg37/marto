#include "kulisch_dim.hpp"


template<int N>
ap_uint<N> acc_IEEE_rounding(
		KulischAcc<N> acc
){
	#pragma HLS INLINE

	ap_uint<1> r_s;
	ap_uint<FPDim<N>::WE> r_e;
	ap_uint<FPDim<N>::WF> r_m;

	r_s = acc[FPDim<N>::ACC_SIZE-1];

	ap_uint<(Static_Val<FPDim<N>::ACC_SIZE>::_log2 +  FPDim<N>::ACC_SIZE)> generic_lzoc = generic_lzoc_shifter<FPDim<N>::ACC_SIZE>(acc, r_s);


	ap_uint<Static_Val<FPDim<N>::ACC_SIZE>::_log2> lzoc = generic_lzoc.range((Static_Val<FPDim<N>::ACC_SIZE>::_log2 +  FPDim<N>::ACC_SIZE)-1, FPDim<N>::ACC_SIZE);
	ap_uint<FPDim<N>::ACC_SIZE> shifted = generic_lzoc.range(FPDim<N>::ACC_SIZE-1, 0);
	

	ap_uint<FPDim<N>::WF> r_m_signed;

	ap_uint<(FPDim<N>::ACC_SIZE - FPDim<N>::WF)> sticky_bits = shifted.range(FPDim<N>::ACC_SIZE-1-FPDim<N>::WF +1-1-1+1,
								FPDim<N>::ACC_SIZE-1-FPDim<N>::WF +1-1-1 - (FPDim<N>::ACC_SIZE - FPDim<N>::WF) +1+1
		                        );


	ap_uint<1> sticky_tmp = sticky_bits.or_reduce();

    ap_uint<1> guard1 = shifted[FPDim<N>::ACC_SIZE-1 -FPDim<N>::WF +1-1-1];
    ap_uint<1> guard2 = shifted[FPDim<N>::ACC_SIZE-1 -FPDim<N>::WF +1-1];

    ap_uint<1> guard, sticky;

	if(lzoc>(FPDim<N>::SUBNORMAL_LIMIT+1))	{
		r_e = 0;
		r_m_signed = shifted.range(FPDim<N>::ACC_SIZE-1, FPDim<N>::ACC_SIZE-1 -FPDim<N>::WF +1) >> (lzoc-(FPDim<N>::SUBNORMAL_LIMIT)-1-1);
		guard=guard2;
		sticky=guard1 or sticky_tmp;

	}
	else{
		r_m_signed = shifted.range(FPDim<N>::ACC_SIZE-1-1, FPDim<N>::ACC_SIZE-1 -FPDim<N>::WF +1-1);
		r_e = FPDim<N>::BIAS-(lzoc)+2;
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

	static constexpr int ceil_log2_acc_size = Static_Val<FPDim<N>::ACC_SIZE>::_log2; 
	static constexpr int pow2_ceil_log2_acc_size = (1<<ceil_log2_acc_size); 
	static constexpr int added_bits = pow2_ceil_log2_acc_size - FPDim<N>::ACC_SIZE;
	ap_uint<1> r_s;
	ap_uint<FPDim<N>::WE> r_e;
	ap_uint<FPDim<N>::WF> r_m;


	ap_int<pow2_ceil_log2_acc_size> ext_acc = (ap_int<FPDim<N>::ACC_SIZE>)acc;

	r_s = acc[FPDim<N>::ACC_SIZE];

	ap_uint<ceil_log2_acc_size + pow2_ceil_log2_acc_size> lzoc_shift = lzoc_shifter<ceil_log2_acc_size>(ext_acc, 0);

	ap_uint<ceil_log2_acc_size> lzoc = lzoc_shift.range(ceil_log2_acc_size + pow2_ceil_log2_acc_size-1, pow2_ceil_log2_acc_size);
	ap_uint<pow2_ceil_log2_acc_size> shifted = lzoc_shift.range(pow2_ceil_log2_acc_size-1, 0);
	
	ap_uint<FPDim<N>::WF> r_m_signed;

	ap_uint<(FPDim<N>::ACC_SIZE - FPDim<N>::WF)> sticky_bits = shifted.range(pow2_ceil_log2_acc_size-1-FPDim<N>::WF +1-1-1,
								pow2_ceil_log2_acc_size-1-FPDim<N>::WF +1-1-1 - (FPDim<N>::ACC_SIZE - FPDim<N>::WF) +1
		                        );

	ap_uint<1> sticky_tmp = sticky_bits.or_reduce();

    ap_uint<1> guard1 = shifted[pow2_ceil_log2_acc_size-1 -FPDim<N>::WF +1-1-1];
    ap_uint<1> guard2 = shifted[pow2_ceil_log2_acc_size-1 -FPDim<N>::WF +1-1];

    ap_uint<1> guard, sticky;

	if(lzoc>(FPDim<N>::SUBNORMAL_LIMIT+added_bits+1))	{
		r_e = 0;
		r_m_signed = shifted.range(pow2_ceil_log2_acc_size-1, pow2_ceil_log2_acc_size-1 -FPDim<N>::WF +1) >> (lzoc-(FPDim<N>::SUBNORMAL_LIMIT+added_bits)-1-1);
		guard=guard2;
		sticky=guard1 or sticky_tmp;
	}
	else{
		r_m_signed = shifted.range(pow2_ceil_log2_acc_size-1-1, pow2_ceil_log2_acc_size-1 -FPDim<N>::WF +1-1);
		r_e = FPDim<N>::BIAS-(lzoc-added_bits)+2;
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
	// ap_uint<1> overflow = (r_m_rounded==0);

	r_m = r_m_rounded_cut;
	
	r_e+=(r_m_rounded[FPDim<N>::WF]);
	ap_uint<1+FPDim<N>::WE> signed_exp = r_s.concat(r_e);


	return signed_exp.concat(r_m);

}