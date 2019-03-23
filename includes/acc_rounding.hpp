#include "kulisch_dim.hpp"


template<int N>
ap_uint<N> acc_IEEE_rounding(
		KulischAcc<N> acc
){
	#pragma HLS INLINE

	static constexpr int ceil_log2_acc_size = Static_Val<FPDim<N>::ACC_SIZE>::_log2; 
	static constexpr int added_bits = (1<<ceil_log2_acc_size) - FPDim<N>::ACC_SIZE;
	ap_uint<1> r_s;
	ap_uint<FPDim<N>::WE> r_e;
	ap_uint<FPDim<N>::WF> r_m;


	ap_int<(1<<ceil_log2_acc_size)> ext_acc = (ap_int<FPDim<N>::ACC_SIZE>)acc;

	r_s = acc[FPDim<N>::ACC_SIZE-1];

	ap_uint<ceil_log2_acc_size + (1<<ceil_log2_acc_size)> lzoc_shift = lzoc_shifter<ceil_log2_acc_size>(ext_acc, r_s);

	ap_uint<ceil_log2_acc_size> lzoc = lzoc_shift.range(ceil_log2_acc_size + (1<<ceil_log2_acc_size)-1, (1<<ceil_log2_acc_size));
	ap_uint<(1<<ceil_log2_acc_size)> shifted = lzoc_shift.range((1<<ceil_log2_acc_size)-1, 0);
	
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

	ap_uint<FPDim<N>::WF+1> r_m_rounded;
	if((guard and not(sticky) and not(r_m_signed[0])) or (guard and sticky)){
		r_m_rounded = r_m_signed+1;
	}
	else{
		r_m_rounded = r_m_signed;	
	}

	ap_uint<FPDim<N>::WF+1> r_m_rounded_cut = r_m_rounded.range(FPDim<N>::WF, 0);
	ap_uint<1> overflow = (r_m_rounded==0) and r_s;
	if(r_s){
		r_m = ~r_m_rounded_cut+1;
	}
	else{
		r_m = r_m_rounded_cut;
	}
	r_e+=(r_m_rounded[FPDim<N>::WF] and not(r_s)) or overflow;
	ap_uint<1+FPDim<N>::WE> signed_exp = r_s.concat(r_e);


	return signed_exp.concat(r_m);

}