#ifndef K1_HPP
#define K1_HPP
#include "kulisch_dim.hpp"

#ifdef K1_ACC_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
inline KulischAcc<WE, WF, Wrapper> add_2CK1(
		KulischAcc<WE, WF, Wrapper> const & acc,
		FPProd<WE, WF, Wrapper> const & prod
){
	constexpr unsigned int PSWidth = IEEEDim<WE, WF>::WFF_Prod;
	constexpr unsigned int ACCWidth = IEEEDim<WE, WF>::ACC_SIZE;
	auto psb = prod.getSignBit();
	auto significand = prod.getSignificand();
	auto inv = significand.invert();
	auto neg = inv.modularAdd(Wrapper<PSWidth, false>{1});
	auto to_add = Wrapper<PSWidth, false>::mux(psb, neg, significand);
	auto sign_ext = Wrapper<ACCWidth - PSWidth, false>::generateSequence(psb);
	auto shifted = sign_ext.concatenate(to_add) << prod.getExp();
	auto r_acc = static_cast<Wrapper<ACCWidth, false> const &>(acc).modularAdd(shifted);
#ifdef K1_ACC_DEBUG
	cerr << "=== K1_ACC ===" << endl;
	cerr << "psb: " << to_string(psb) << endl;
	cerr << "significand: " << to_string(significand) << endl;
	cerr << "inv: " << to_string(inv) << endl;
	cerr << "neg: " << to_string(neg) << endl;
	cerr << "to_add: " << to_string(to_add) << endl;
	cerr << "sign_ext: " << to_string(sign_ext) << endl;
	cerr << "shifted: " << to_string(shifted) << endl;
	cerr << "r_acc: " << to_string(r_acc) << endl;
	cerr << "==============" << endl;
#endif
	return r_acc;
}
#endif
