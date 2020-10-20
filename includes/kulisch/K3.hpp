
#include "kulisch_dim.hpp"
using namespace std;
#include "bitvector/lzoc_shifter.hpp"
#include "bitvector/shifter_marto.hpp"

using hint::Static_Val;
using hint::Static_Ceil_Div;

template <unsigned int stageIndex, unsigned int spread, unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
	Wrapper<WE+1 +1 - Static_Val<bankSize>::_log2 +1, false> const & stageSelect,
	Wrapper<1, false> const & inputSign,
	Wrapper<(1<<Static_Val<spread*bankSize>::_log2), false> const &,
	typename enable_if<(spread < 1)>::type* = 0
)
{
	constexpr unsigned int MantSpread = getMantSpread<WE, WF, bankSize>();

	auto cond = (stageSelect + MantSpread) < stageIndex;
	auto ret = Wrapper<bankSize, false>::generateSequence(cond.invert() & inputSign);
	return ret;
}

template <unsigned int stageIndex, unsigned int spread, unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
	Wrapper<WE+1 +1 - Static_Val<bankSize>::_log2 +1, false> const & stageSelect,
	Wrapper<1, false> const & inputSign,
	Wrapper<(1<<Static_Val<spread*bankSize>::_log2), false> const & shiftedSignificand,
	typename enable_if<(spread >= 1)>::type* = 0
	)
{
	constexpr unsigned int select_size = WE+ 1 + 1 - Static_Val<bankSize>::_log2 + 1;
	constexpr unsigned int ident_x = stageIndex - spread + 1;
	constexpr unsigned int up = spread * bankSize - 1;
	constexpr unsigned int down = (spread - 1)*bankSize;

	Wrapper<select_size, false> activationValue{ident_x};

	auto cond = (stageSelect == activationValue);

	auto res = Wrapper<bankSize, false>::mux(
				cond,
				shiftedSignificand.template slice<up, down>(),
				add_2CK3_to_add_Rec<stageIndex, spread-1>(stageSelect, inputSign, shiftedSignificand)
		);

	return res;
}


template <unsigned int stageIndex, unsigned int spread, unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline Wrapper<bankSize, false> add_2CK3_to_add_(
	Wrapper<WE+1 +1 - Static_Val<bankSize>::_log2 +1, false> const & stageSelect,
	Wrapper<1, false> const & inputSign,
	Wrapper<(1<<Static_Val<spread*bankSize>::_log2), false> const & shiftedSignificand
		) {
	return add_2CK3_to_add_Rec<stageIndex, spread>(stageSelect, inputSign, shiftedSignificand);
}

template<unsigned int stageIndex, unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline Wrapper<bankSize+1, false> add_2CK3_acc_stage(
					acc_2CK3<WE, WF, bankSize, Wrapper> const & acc,
					Wrapper<WE+1 +1 - Static_Val<bankSize>::_log2 +1, false> const & stageSelect,
					Wrapper<(1<<Static_Val<getMantSpread<WE, WF, bankSize>()*bankSize>::_log2), false> const & shiftedSignificand,
					Wrapper<1, false> const & sign,
					Wrapper<1, false> const & carry0 = {0}
		)
{
	auto bank = acc.template getBank<stageIndex>();
	auto accCarry = (stageIndex==0) ? carry0 : acc.template getCarry<stageIndex-1>();

	auto toAdd = add_2CK3_to_add_<WE, WF, bankSize, getMantSpread<WE, WF, bankSize>(), stageIndex>(stageSelect, sign, shiftedSignificand);
	return bank.addWithCarry(toAdd, accCarry);
}

template<unsigned int stage, unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline void _perform_addition_all_stages(
			acc_2CK3<WE, WF, bankSize, Wrapper>const & acc,
			Wrapper<WE+1 +1 - Static_Val<bankSize>::_log2 +1, false> const & stageSelect,
			Wrapper<(1<<Static_Val<getMantSpread<WE, WF, bankSize>()*bankSize>::_log2), false> const & shiftedSignificand,
			Wrapper<1, false> const & sign,
			acc_2CK3<WE, WF, bankSize, Wrapper>& resAcc,
			typename enable_if<(stage == 0)>::type* = 0
		)
{
	auto stageResult = add_2CK3_acc_stage<stage>(
				acc,
				stageSelect,
				shiftedSignificand,
				sign
	);
	resAcc.template setCarry<stage>(stageResult.template get<bankSize>());
	resAcc.template setBank<stage>(stageResult.template slice<bankSize-1,0>());
}

template<unsigned int stage, unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline void _perform_addition_all_stages(
			acc_2CK3<WE, WF, bankSize, Wrapper>const & acc,
			Wrapper<WE+1 +1 - Static_Val<bankSize>::_log2 +1, false> const & stageSelect,
			Wrapper<(1<<Static_Val<getMantSpread<WE, WF, bankSize>()*bankSize>::_log2), false> const & shiftedSignificand,
			Wrapper<1, false> const & sign,
			acc_2CK3<WE, WF, bankSize, Wrapper>& resAcc,
			typename enable_if<(stage > 0)>::type* = 0
		)
{
	auto stageResult = add_2CK3_acc_stage<stage>(
				acc,
				stageSelect,
				shiftedSignificand,
				sign
	);
	resAcc.template setCarry<stage>(stageResult.template get<bankSize>());
	resAcc.template setBank<stage>(stageResult.template slice<bankSize-1,0>());
	_perform_addition_all_stages<stage-1>(acc, stageSelect, shiftedSignificand, sign, resAcc);
}

template<unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline acc_2CK3<WE, WF, bankSize, Wrapper> add_2CK3(
		acc_2CK3<WE, WF, bankSize, Wrapper> const & acc,
		FPProd<WE, WF, Wrapper> const & prod
	)
{
	constexpr unsigned int spread = Static_Ceil_Div<2*WF+2,bankSize>::val+1;
	constexpr unsigned int bits_for_shift = Static_Val<spread*bankSize>::_log2;
	constexpr unsigned int nb_stages = getNbStages<WE, WF, bankSize>();

	auto prodExp = prod.getExp().template leftpad<IEEEProdDim<WE, WF>::WProdExp + 1>();
	auto prodSign = prod.getSignBit();

	auto inputSignificand = prod.getSignificand();
	auto inputSignificandComplemented = Wrapper<IEEEProdDim<WE, WF>::WProdFrac + 1, false>::mux(
				prod.getSignBit(),
				(inputSignificand.invert() + Wrapper<IEEEProdDim<WE, WF>::WProdFrac, false>{1}),
				 prodExp.template leftpad<IEEEProdDim<WE, WF>::WProdFrac + 1>()
		);

	auto inputSignificandWithSign = prodSign.concat(inputSignificandComplemented);
	auto shiftValue = prodExp.range(Static_Val<bankSize>::_log2-1,0);
	auto ext = inputSignificandWithSign.template leftpad<(1<<bits_for_shift)>();

	auto shiftedInput = shifter<bits_for_shift>(ext,shiftValue,0);

	auto stageSelect = prodExp.template slice<WE + 1, Static_Val<bankSize>::_log2>();
	acc_2CK3<WE, WF, bankSize, Wrapper> fullAcc{};

	_perform_addition_all_stages<nb_stages-1>(acc, stageSelect, shiftedInput, prodSign, fullAcc);

	return fullAcc;
}

template<unsigned int nb_repeat, unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline acc_2CK3<WE, WF, bankSize, Wrapper>  _propagate_carries_2CK3_Req(
		acc_2CK3<WE, WF, bankSize, Wrapper> const & acc,
		typename enable_if<(nb_repeat > 0)>::type* = 0
)
{
		constexpr unsigned int nb_stages = getNbStages<WE, WF, bankSize>();
		acc_2CK3<WE, WF, bankSize, Wrapper> ret;
		_perform_addition_all_stages<nb_stages - 1>(acc, {{0}}, {{0}}, ret);
		return _propagate_carries_2CK3_Req<nb_repeat - 1>(ret);
}

template<unsigned int nb_repeat, unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline acc_2CK3<WE, WF, bankSize, Wrapper>  _propagate_carries_2CK3_Req(
		acc_2CK3<WE, WF, bankSize, Wrapper> const & acc,
		typename enable_if<(nb_repeat == 0)>::type* = 0
)
{
		constexpr unsigned int nb_stages = getNbStages<WE, WF, bankSize>();
		acc_2CK3<WE, WF, bankSize, Wrapper> ret;
		_perform_addition_all_stages<nb_stages - 1>(acc, {{0}}, {{0}}, ret);
		return ret;
}

template <unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline KulischAcc<WE, WF, Wrapper> propagate_carries_2CK3(acc_2CK3<WE, WF, bankSize, Wrapper> const & acc)
{
	constexpr unsigned int nb_stages = getNbStages<WE, WF, bankSize>();
	auto ret = _propagate_carries_2CK3_Req<nb_stages>(acc);
	return ret.getAcc();
}
