#pragma once
#include <cstdio>
#include <utility>

#include "posit_dim.hpp"
#include "primitives/lzoc_shifter.hpp"

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY>
inline Quire<N, WES, Wrapper, NB_CARRY> add_sub_quire(
		Quire<N, WES, Wrapper, NB_CARRY> quire,
		PositProd<N, WES, Wrapper> input,
		Wrapper<1, false> isSub
){
	constexpr int LOG2_EXT_SUM_SIZE = hint::Static_Val<quire.Size-1 +input.SignificandSize>::_log2;
	constexpr int QUIRE_SIZE = quire.Size;
	constexpr int PROD_SIGNIFICAND_SIZE = PositDim<N, WES>::ProdSignificandSize;

	auto inputSignificand = input.getSignedSignificand();
	auto sign = input.getSignBit();
	auto replicated_sign = Wrapper<PROD_SIGNIFICAND_SIZE + 1, false>::generateSequence(sign);
	auto complementedInputIfIsSub = replicated_sign ^ inputSignificand;

	auto shiftValue = input.getExp();
	constexpr unsigned int SHIFT_SIZE  = 1<<LOG2_EXT_SUM_SIZE;
	auto ext = complementedInputIfIsSub.template leftpad<SHIFT_SIZE>();

	auto shiftedInput = hint::shifter<false, SHIFT_SIZE, input.ExpSize>(ext, shiftValue, isSub);
	auto shiftedInputShrinked = shiftedInput.template slice<
				QUIRE_SIZE + PROD_SIGNIFICAND_SIZE-2,
				PROD_SIGNIFICAND_SIZE
		>();

	auto quireWithoutSignAndNARBit = quire.getQuireWithoutNaR();
	auto sumResult = shiftedInputShrinked.addWithCarry(
			quireWithoutSignAndNARBit,
			isSub
		).template slice<quire.Size-2, 0>();

	auto resultIsNaR = quire.getIsNaR() | input.getIsNaR();

	return Quire<N, WES, Wrapper, NB_CARRY>{resultIsNaR.concatenate(sumResult)};
}


template<unsigned int bankSize>
constexpr unsigned int getIndex(int index, bool isUpper)
{
	return bankSize*index - isUpper;
}

template <unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int bankSize, unsigned int spread, unsigned int loop_idx>
inline Wrapper<bankSize, false> getToAddRec(
	Wrapper<PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>(), false> stageSelect,
	Wrapper<1, false> inputSign,
	Wrapper<1, false> isSub,
	Wrapper<getExtShiftSize<N, WES, bankSize>(), false>,
	typename enable_if<(spread == 0)>::type* = 0
)
{
	constexpr unsigned int SS_WIDTH = PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>();
	constexpr unsigned int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositDim<N, WES>::ProdSignificandSize, bankSize>::val;
	constexpr unsigned int BOUND = loop_idx + BANKS_FOR_USELESS_BITS - getMantSpread<N, WES, bankSize>();
	Wrapper<SS_WIDTH + 1, false> bound_wrapper{{BOUND}};
	Wrapper<SS_WIDTH, false> loop_idx_wrap {{loop_idx}};
	auto test_wrap = (stageSelect.template leftpad<SS_WIDTH + 1>() <  bound_wrapper);
	auto loop_idx_below_wrap = loop_idx_wrap < stageSelect;
	Wrapper<1, false> fill_bit = (test_wrap & (inputSign xor isSub)) | (loop_idx_below_wrap & isSub);
	return Wrapper<bankSize, false>::generateSequence(fill_bit);
}

template <unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int bankSize, unsigned int spread, unsigned int loop_idx>
inline Wrapper<bankSize, false> getToAddRec(
	Wrapper<PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>(), false> stageSelect,
	Wrapper<1, false> inputSign,
	Wrapper<1, false> isSub,
	Wrapper<getExtShiftSize<N, WES, bankSize>(), false> shiftedSignificand,
	typename enable_if<(spread > 0)>::type* = 0
)
{
	constexpr unsigned int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositDim<N, WES>::ProdSignificandSize, bankSize>::val;
	// We want to check if loop_idx == (stageSelect+spread-1-BANKS_FOR_USELESS_BITS)
	// Transformed in stageSelect == loop_idx - spread + 1 + BANK_FOR_USELESS_BITS
	constexpr unsigned int test_val = loop_idx + 1 + BANKS_FOR_USELESS_BITS - spread;
	Wrapper<PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>(), false> test_val_wrap{{test_val}};
	auto check = (test_val_wrap == stageSelect);
	auto ret = Wrapper<bankSize, false>::mux(
			check,
			shiftedSignificand.template slice<spread*bankSize-1,(spread-1)*bankSize>(),
			getToAddRec<N, WES, Wrapper, NB_CARRY, bankSize, (spread-1), loop_idx>(stageSelect, inputSign, isSub, shiftedSignificand)
		);
	return ret;
}


template <unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int bankSize, unsigned int spread, unsigned int loop_idx>
inline Wrapper<bankSize, false> getToAdd(
	Wrapper<PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>(), false> stageSelect,
	Wrapper<1, false> inputSign,
	Wrapper<1, false> isSub,
	Wrapper<getExtShiftSize<N, WES, bankSize>(), false> shiftedSignificand
		){
	return getToAddRec<N, WES, Wrapper, NB_CARRY, bankSize, spread, loop_idx>(stageSelect, inputSign, isSub, shiftedSignificand);
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int bankSize, unsigned int loop_idx>
inline Wrapper<bankSize+1, false> add_sub_quire_stage(SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> quire,
										Wrapper<PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>(), false> stageSelect,
										Wrapper<1, false> inputSign,
										Wrapper<1, false> isSub,
										Wrapper<getExtShiftSize<N, WES, bankSize>(), false> shiftedSignificand,
										typename enable_if<(loop_idx > 0)>::type* = 0
	)
{
	auto quireBank = quire.template getBank<loop_idx>();
	Wrapper<1, false> quireCarry = quire.template getCarry<loop_idx - 1>();

	auto toAdd = getToAdd<N, WES, Wrapper, NB_CARRY, bankSize, getMantSpread<N, WES, bankSize>(), loop_idx>(stageSelect, inputSign, isSub, shiftedSignificand);

	auto sum = quireBank.addWithCarry(toAdd, quireCarry);
	return sum;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int bankSize, unsigned int loop_idx>
inline Wrapper<bankSize+1, false> add_sub_quire_stage(SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> quire,
										Wrapper<PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>(), false> stageSelect,
										Wrapper<1, false> inputSign,
										Wrapper<1, false> isSub,
										Wrapper<getExtShiftSize<N, WES, bankSize>(), false> shiftedSignificand,
										typename enable_if<(loop_idx == 0)>::type* = 0
	)
{
	auto quireBank = quire.template getBank<loop_idx>();
	Wrapper<1, false> quireCarry = isSub;

	auto toAdd = getToAdd<N, WES, Wrapper, NB_CARRY, bankSize, getMantSpread<N, WES, bankSize>(), loop_idx>(stageSelect, inputSign, isSub, shiftedSignificand);

	auto sum = quireBank.addWithCarry(toAdd, quireCarry);
	return sum;
}

/**
 * TODO :
 *
 *  - Convert add_sub_quire_stage
 *  - Loop flatening using template unrolling / concatenation
 */

template<unsigned int loop_idx, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
using _partial_seg_quire_store = std::pair<Wrapper<(loop_idx + 1)*bankSize, false>, Wrapper<loop_idx + 1, false>>;

template<unsigned int NB_CARRY, unsigned int bankSize, unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int loop_idx>
inline _partial_seg_quire_store<loop_idx, bankSize, Wrapper> _seq_quire_partial_add(
			SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> quire,
			Wrapper<PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>(), false> stage_select,
			Wrapper<1, false> sign,
			Wrapper<1, false> isSub,
			Wrapper<getExtShiftSize<N, WES, bankSize>(), false> input,
			typename enable_if<(loop_idx > 0)>::type* = 0
		)
{
	auto stageResult = add_sub_quire_stage<N, WES, Wrapper, NB_CARRY, bankSize, loop_idx>(
				quire,
				stage_select,
				sign,
				isSub,
				input
		);
	auto lower_stages = _seq_quire_partial_add<NB_CARRY, bankSize, N, WES, Wrapper, loop_idx - 1>(quire, stage_select, sign, isSub, input);
	_partial_seg_quire_store<loop_idx, bankSize, Wrapper> ret = make_pair(
			lower_stages.first.concatenate(stageResult.template slice<bankSize - 1, 0>()),
			lower_stages.second.concatenate(stageResult.template get<bankSize>())
		);
	return ret;
}

template<unsigned int NB_CARRY, unsigned int bankSize, unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int loop_idx>
inline _partial_seg_quire_store<loop_idx, bankSize, Wrapper> _seq_quire_partial_add(
			SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> quire,
			Wrapper<PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>(), false> stage_select,
			Wrapper<1, false> sign,
			Wrapper<1, false> isSub,
			Wrapper<getExtShiftSize<N, WES, bankSize>(), false> input,
			typename enable_if<(loop_idx== 0)>::type* = 0
		)
{
	 Wrapper<bankSize+1, false> stageResult = add_sub_quire_stage<N, WES, Wrapper, NB_CARRY, bankSize, 0>(
				quire,
				stage_select,
				sign,
				isSub,
				input
		);
	_partial_seg_quire_store<loop_idx, bankSize, Wrapper> ret{stageResult.template slice<bankSize - 1, 0>(), stageResult.template get<bankSize>()};
	return ret;
}


template<unsigned int NB_CARRY, unsigned int bankSize, unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<QuireDim<N, WES, NB_CARRY>::Size + getNbStages<N, WES, NB_CARRY, bankSize>() - 1, false>
_perform_seg_add(
		SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> quire,
		Wrapper<PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>(), false> stage_select,
		Wrapper<1, false> sign,
		Wrapper<1, false> isSub,
		Wrapper<getExtShiftSize<N, WES, bankSize>(), false> input
		)
{
	constexpr unsigned int start_idx = getNbStages<N, WES, NB_CARRY, bankSize>() - 1;
	auto val = _seq_quire_partial_add<NB_CARRY, bankSize, N, WES, Wrapper, start_idx>(
			quire,
			stage_select,
			sign,
			isSub,
			input
		);
	return val.first.concatenate(val.second);
}

template<unsigned int NB_CARRY, unsigned int bankSize, unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> segmented_add_sub_quire(
		SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> quire,
		PositProd<N, WES, Wrapper> input,
		Wrapper<1, false> isSub)
{
	//static constexpr int SHIFTED_SIGNIFICAND_SIZE = PositProd<N, WES>::SignificandSize+1 + (1<<getShiftSize<bankSize>());
	constexpr unsigned int PROD_EXP_SIZE = PositDim<N, WES>::ProdExpSize;
	constexpr unsigned int PROD_SIGIFICAND_SIZE = PositDim<N, WES>::ProdSignificandSize;
	// constexpr unsigned int BANKS_FOR_USELESS_BITS = hint::Static_Ceil_Div<PROD_SIGIFICAND_SIZE, bankSize>::val;
	// constexpr unsigned int padding = bankSize*BANKS_FOR_USELESS_BITS - PROD_SIGIFICAND_SIZE;
	// constexpr unsigned int ENCODING_BITS_FOR_USELESS_BANKS = hint::Static_Val<BANKS_FOR_USELESS_BITS>::_log2;
	constexpr unsigned int LOG2_SHIFT_SIZE = Static_Val<getExtShiftSize<N, WES, bankSize>()>::_log2;
	constexpr unsigned int NB_STAGES = getNbStages<N, WES, NB_CARRY, bankSize>();
	constexpr unsigned int PROD_EXP_BANK_SIZE_MOD = PROD_EXP_SIZE % bankSize;
	constexpr unsigned int SHIFT_SIZE = getShiftSize<bankSize>();
	constexpr unsigned int EXT_SHIFT_SIZE = getExtShiftSize<N, WES, bankSize>();
	constexpr unsigned int STAGE_SELECT_SIZE = PROD_EXP_SIZE - SHIFT_SIZE;

	//Wrapper<PROD_EXP_SIZE+ENCODING_BITS_FOR_USELESS_BANKS, false> prodExp = input.getExp()+padding;
	auto prod_exp = input.getExp();
	auto input_significand = input.getSignedSignificand();
	auto sign = input_significand.template get<PROD_SIGIFICAND_SIZE>();
	auto extended_sign = Wrapper<PROD_SIGIFICAND_SIZE + 1, false>::generateSequence(sign);

	auto complementedInputIfIsSub = extended_sign ^ input_significand;
	auto last_bit_exp_extended = (prod_exp.template slice<SHIFT_SIZE-1,0>()).template leftpad<SHIFT_SIZE+1>();

	Wrapper<SHIFT_SIZE + 1, false> shift_modulus{PROD_EXP_BANK_SIZE_MOD};

	auto last_bit_exp_ext_sub = last_bit_exp_extended.modularSub(shift_modulus);
	auto exp_borrow = last_bit_exp_ext_sub.template get<SHIFT_SIZE>().template leftpad<STAGE_SELECT_SIZE>();
	auto exp_low_bits = last_bit_exp_ext_sub.template slice<SHIFT_SIZE - 1, 0>();

	auto ext = complementedInputIfIsSub.template leftpad<(1<<LOG2_SHIFT_SIZE)>();
	auto shifted_input = hint::shifter<false, (1<<LOG2_SHIFT_SIZE), SHIFT_SIZE, false>(ext, exp_low_bits, isSub);
	auto shifted_input_shrinked = shifted_input.template slice<EXT_SHIFT_SIZE-1, 0>();

	auto stage_select = prod_exp.template slice<PROD_EXP_SIZE-1, SHIFT_SIZE>().modularSub(exp_borrow);

	Wrapper<1, false> resultIsNaR = input.getIsNaR().bitwise_or(quire.getIsNaR());
	auto fullQuireWoNaR = _perform_seg_add(quire, stage_select, sign, isSub, shifted_input_shrinked);
	return SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize>{resultIsNaR.concatenate(fullQuireWoNaR)};
}

template<unsigned int NB_CARRY, unsigned int bankSize, unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int loop_idx>
inline SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize>
_propag_carry_seg_req(
		SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> quire,
		typename enable_if<(loop_idx > 0)>::type* = 0
		)
{
	SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> res {
			quire.getIsNaR().concatenate(
				_perform_seg_add(quire, {0}, {0}, {0}, {0})
			)
		};
	return _propag_carry_seg_req<NB_CARRY, bankSize, N, WES, Wrapper, loop_idx - 1>(res);
}

template<unsigned int NB_CARRY, unsigned int bankSize, unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int loop_idx>
inline SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize>
_propag_carry_seg_req(
		SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> quire,
		typename enable_if<(loop_idx == 0)>::type* = 0
		)
{
	SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> res {
			quire.getIsNaR().concatenate(
				_perform_seg_add(quire, {0}, {0}, {0}, {0})
			)
		};
	return res;
}


template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int bankSize>
inline Quire<N, WES, Wrapper, NB_CARRY> propagateCarries(SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> quire)
{
	constexpr unsigned int NB_STAGES = getNbStages<N, WES, NB_CARRY, bankSize>();
	SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> fullQuire {quire};
	auto res = _propag_carry_seg_req<NB_CARRY, bankSize, N, WES, Wrapper, NB_STAGES - 1>(fullQuire);
	return res.getAsQuireWoCarries();
}
