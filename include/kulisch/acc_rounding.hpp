#ifndef ACC_ROUNDING_HPP
#define ACC_ROUNDING_HPP

#include "primitives/lzoc_shifter.hpp"
#include "primitives/shifter.hpp"

#include "kulisch_dim.hpp"
#include "tools/static_math.hpp"


#ifdef KULISCH_ROUND_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif

using hint::Static_Val;
using hint::LZOC_shift;

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
inline IEEENumber<WE, WF, Wrapper> acc_IEEE_rounding(
		KulischAcc<WE, WF, Wrapper> const & acc
){
	constexpr unsigned int subnormal_first_idx = IEEEDim<WE, WF>::MIN_SN_ACC_OFFSET;
	constexpr unsigned int normal_max_idx = IEEEDim<WE, WF>::MAX_N_ACC_OFFSET;
	constexpr unsigned int overflow_bits = IEEEDim<WE, WF>::ACC_SIZE - (normal_max_idx + 2);

	constexpr unsigned int reduced_acc_size = normal_max_idx - subnormal_first_idx + 2;

	constexpr unsigned int lzoc_in_size = reduced_acc_size;
	constexpr unsigned int lzoc_c2pow = Static_Val<lzoc_in_size>::_2pow;
	constexpr unsigned int lzoc_count = (lzoc_c2pow == lzoc_in_size) ? lzoc_c2pow : lzoc_c2pow - 1;
	constexpr unsigned int lzoc_size = Static_Val<lzoc_count>::_storage;

	//constexpr unsigned int ceil_log2_acc_size_minus_1 = Static_Val<reduced_acc_size>::_rlog2-1;
	//constexpr unsigned int pow2_ceil_log2_acc_size_minus_1 = (1<<ceil_log2_acc_size_minus_1);
	//constexpr unsigned int remaining_bits = reduced_acc_size - (1 << (Static_Val<reduced_acc_size>::_rlog2)-1);
	//constexpr unsigned int ceil_log2_remaining_bits = Static_Val<remaining_bits>::_rlog2;
	//constexpr unsigned int pow2_ceil_log2_remaining_bits = (1<<ceil_log2_remaining_bits);

	// Exp bit position is exponent + 2*Bias so under 1 bias we have only
	// We keep one extra bit for having the real rounding value for subnormals
	auto low_acc = acc.template slice<subnormal_first_idx - 2 , 0>();
	auto mid_acc = acc.template slice<normal_max_idx, subnormal_first_idx - 1>();
	auto high_acc = acc.template slice<IEEEDim<WE, WF>::ACC_SIZE -2, normal_max_idx + 1>();

	auto r_s = acc.getSignBit();

	auto overflow_when_neg = high_acc.bitwise_xor(Wrapper<overflow_bits, false>::generateSequence(r_s)).or_reduction();
	auto sticky_low = low_acc.or_reduction();
	// Important : we concatenate RS to allow good 2's complement computation of subnormal negative 2 powers
	auto lzoc_shift = LZOC_shift<lzoc_in_size, lzoc_count>(mid_acc, r_s);
	auto lzoc = lzoc_shift.lzoc;
	auto shifted = lzoc_shift.shifted;

	auto sticky_bits = shifted.template slice<reduced_acc_size - 3 - WF, 0>();

	auto sticky_tmp = sticky_bits.or_reduction() | sticky_low;

	auto guard1 = shifted.template get<reduced_acc_size -2 - WF>();
	auto guard2 = shifted.template get<reduced_acc_size -1 - WF>();
	auto guard3 = acc.template get<subnormal_first_idx-1>();

	constexpr unsigned int lzoc_bound = reduced_acc_size - (WF + 1);
	constexpr unsigned int lzoc_decr = lzoc_bound + 1;
	Wrapper<lzoc_size, false> lzoc_bound_wrap {lzoc_bound}, lzoc_decr_wrap{lzoc_decr}, limit_subnormal_wrap{lzoc_in_size - 1};
	auto is_subnormal = lzoc > lzoc_bound_wrap;

	//TODO Revoir size r_m_signed etc. à partir d'ici

	//Subnormal case
	Wrapper<WE, false> sn_exp{0};
	auto sn_shift = lzoc.modularSub(lzoc_decr_wrap);
	auto sn_shift_slice = shifted.template slice<reduced_acc_size-1, reduced_acc_size - WF>();
	auto sn_signed_signif = hint::shifter<true>(sn_shift_slice, sn_shift, r_s);
	auto sn_limit = lzoc == limit_subnormal_wrap;
	auto sn_guard = Wrapper<1, false>::mux(sn_limit, guard3, guard2);
	auto sn_sticky = Wrapper<1, false>::mux(sn_limit, sticky_tmp, guard1 | guard3 | sticky_tmp);

	//Normal case
	constexpr unsigned int bias = IEEEDim<WE, WF>::MAX_NORMAL_BIASED_EXP;
	Wrapper<lzoc_size, false> bias_wrap{bias};
	auto n_signed_signif = shifted.template slice<reduced_acc_size - 2, reduced_acc_size - WF - 1>();
	auto n_exp = bias_wrap.modularSub(lzoc).template slice<WE-1, 0>();
	auto n_guard = guard1;
	auto n_sticky = sticky_tmp | guard3;

	// Affectation
	auto r_e = Wrapper<WE, false>::mux(is_subnormal, sn_exp, n_exp);
	auto r_m_signed = Wrapper<WF, false>::mux(is_subnormal, sn_signed_signif, n_signed_signif);
	auto guard = Wrapper<1, false>::mux(is_subnormal, sn_guard, n_guard);
	auto sticky = Wrapper<1, false>::mux(is_subnormal, sn_sticky, n_sticky);
	auto low_bit = r_m_signed.template get<0>() ^ r_s;

	auto roundBit = (guard & ((sticky.invert() & low_bit) | (sticky & r_s.invert()))) | (r_s & guard.invert());

	auto r_m_rounded = r_m_signed.addWithCarry(roundBit.template leftpad<WF>(), {{0}});
	auto r_m = Wrapper<WF, false>::mux(r_s,
									   r_m_signed.invert(),
									   r_m_signed);

	auto pre_round = r_s.concatenate(r_e).concatenate(r_m);
	auto res = pre_round.modularAdd(roundBit.template leftpad<WE+WF+1>());
#ifdef KULISCH_ROUND_DEBUG
	cerr << "=== KULISCH_ROUND ===" << endl;
	cerr << "acc: " << to_string(acc.downcast()) << endl;
	cerr << "low_acc: " << to_string(low_acc) << endl;
	cerr << "mid_acc: " << to_string(mid_acc) << endl;
	cerr << "high_acc: " << to_string(high_acc) << endl;
	cerr << "r_s: " << to_string(r_s) << endl;
	cerr << "overflow_when_neg: " << to_string(overflow_when_neg) << endl;
	cerr << "sticky_low: " << to_string(sticky_low) << endl;
	cerr << "lzoc: " << to_string(lzoc) << endl;
	cerr << "shifted: " << to_string(shifted) << endl;
	cerr << "sticky_bits: " << to_string(sticky_bits) << endl;
	cerr << "sticky_tmp: " << to_string(sticky_tmp) << endl;
	cerr << "guard1: " << to_string(guard1) << endl;
	cerr << "guard2: " << to_string(guard2) << endl;
	cerr << "guard3: " << to_string(guard3) << endl;
	cerr << "is_subnormal: " << to_string(is_subnormal) << endl;
	cerr << "sn_shift: " << to_string(sn_shift) << endl;
	cerr << "sn_shift_slice: " << to_string(sn_shift_slice) << endl;
	cerr << "sn_signed_signif: " << to_string(sn_signed_signif) << endl;
	cerr << "sn_limit: " << to_string(sn_limit) << endl;
	cerr << "sn_guard: " << to_string(sn_guard) << endl;
	cerr << "sn_sticky: " << to_string(sn_sticky) << endl;
	cerr << "n_signed_signif: " << to_string(n_signed_signif) << endl;
	cerr << "n_exp: " << to_string(n_exp) << endl;
	cerr << "n_guard: " << to_string(n_guard) << endl;
	cerr << "n_sticky: " << to_string(n_sticky) << endl;
	cerr << "r_e: " << to_string(r_e) << endl;
	cerr << "r_m_signed: " << to_string(r_m_signed) << endl;
	cerr << "guard: " << to_string(guard) << endl;
	cerr << "sticky: " << to_string(sticky) << endl;
	cerr << "low_bit: " << to_string(low_bit) << endl;
	cerr << "roundBit: " << to_string(roundBit) << endl;
	cerr << "r_m_rounded: " << to_string(r_m_rounded) << endl;
	cerr << "r_m: " << to_string(r_m) << endl;
	cerr << "pre_round: " << to_string(pre_round) << endl;
	cerr << "res: " << to_string(res) << endl;
	cerr << "====================" << endl;
#endif
	return res;
}
#endif

