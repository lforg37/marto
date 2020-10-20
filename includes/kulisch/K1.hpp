#include "kulisch_dim.hpp"

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
inline KulischAcc<WE, WF, Wrapper> add_2CK1(
		KulischAcc<WE, WF, Wrapper> const & acc,
		FPProd<WE, WF, Wrapper> const & prod
){
	constexpr unsigned int PSWidth = IEEEProdDim<WE, WF>::WProdFrac;
	constexpr unsigned int ACCWidth = FPDim<WE, WF>::ACC_SIZE;
	auto ext_significand = prod.getSignificand();
	auto psb = prod.getSignBit();
	auto significand = prod.getSignificand();
	auto inv = significand.invert();
	auto neg = inv + Wrapper<PSWidth, false>{1};
	auto to_add = Wrapper<PSWidth + 1, false>::mux(psb, neg, significand.template leftpad<PSWidth + 1>());
	auto sign_ext = Wrapper<ACCWidth - PSWidth - 1, false>::generateSequence(to_add.template get<PSWidth>());
	auto shifted = sign_ext.concatenate(to_add) << prod.getExp();
	auto r_acc = shifted + acc;
	return r_acc;
}
