#pragma once
#include <cstdio>
#include <utility>

#include "posit_dim.hpp"
#include "primitives/lzoc_shifter.hpp"

#ifdef POSIT_QUIRE_ADD_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif

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
	auto replicated_sign = Wrapper<PROD_SIGNIFICAND_SIZE + 1, false>::generateSequence(isSub);
	auto complementedInputIfIsSub = replicated_sign ^ inputSignificand;

	auto shiftValue = input.getExp().modularAdd({{2*PositDim<N, WES>::EMax + 1}});
	constexpr unsigned int SHIFT_SIZE  = 1<<LOG2_EXT_SUM_SIZE;
	auto ext = complementedInputIfIsSub.as_signed().template leftpad<SHIFT_SIZE>().as_unsigned();

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

#ifdef POSIT_QUIRE_ADD_DEBUG
	cerr << "=== Add Sub Quire ===" << endl;
	cerr << "inputSignificand: " << to_string(inputSignificand) << endl;
	cerr << "sign: " << to_string(sign) << endl;
	cerr << "replicated_sign: " << to_string(replicated_sign) << endl;
	cerr << "complementedInputIfIsSub: " << to_string(complementedInputIfIsSub) << endl;
	cerr << "shiftValue: " << to_string(shiftValue) << endl;
	cerr << "ext: " << to_string(ext) << endl;
	cerr << "shiftedInput: " << to_string(shiftedInput) << endl;
	cerr << "shiftedInputShrinked: " << to_string(shiftedInputShrinked) << endl;
	cerr << "quireWithoutSignAndNARBit: " << to_string(quireWithoutSignAndNARBit) << endl;
	cerr << "sumResult: " << to_string(sumResult) << endl;
	cerr << "resultIsNaR: " << to_string(resultIsNaR) << endl;
	cerr << "==============================" << endl;
#endif

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
	constexpr int TMP_BOUND = loop_idx + BANKS_FOR_USELESS_BITS - getQuireSpread<N, WES, bankSize>();
	constexpr unsigned int BOUND = (TMP_BOUND < 0) ? 0 : TMP_BOUND;
	Wrapper<SS_WIDTH + 1, false> bound_wrapper{{BOUND}};
	Wrapper<SS_WIDTH, false> loop_idx_wrap {{loop_idx}};
	auto test_wrap = (stageSelect.template leftpad<SS_WIDTH + 1>() <=  bound_wrapper);
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

	auto toAdd = getToAdd<N, WES, Wrapper, NB_CARRY, bankSize, getQuireSpread<N, WES, bankSize>(), loop_idx>(stageSelect, inputSign, isSub, shiftedSignificand);


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

	auto toAdd = getToAdd<N, WES, Wrapper, NB_CARRY, bankSize, getQuireSpread<N, WES, bankSize>(), loop_idx>(stageSelect, inputSign, isSub, shiftedSignificand);

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
	auto stage_result = add_sub_quire_stage<N, WES, Wrapper, NB_CARRY, bankSize, loop_idx>(
				quire,
				stage_select,
				sign,
				isSub,
				input
		);
	auto lower_stages = _seq_quire_partial_add<NB_CARRY, bankSize, N, WES, Wrapper, loop_idx - 1>(quire, stage_select, sign, isSub, input);
	auto cur_stage_carry = stage_result.template get<bankSize>();
	auto cur_stage_bank = stage_result.template slice<bankSize - 1, 0>();
	_partial_seg_quire_store<loop_idx, bankSize, Wrapper> ret = make_pair(
			cur_stage_bank.concatenate(lower_stages.first),
			cur_stage_carry.concatenate(lower_stages.second)
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
	 auto stage_result = add_sub_quire_stage<N, WES, Wrapper, NB_CARRY, bankSize, 0>(
				quire,
				stage_select,
				sign,
				isSub,
				input
		);
	_partial_seg_quire_store<loop_idx, bankSize, Wrapper> ret{stage_result.template slice<bankSize - 1, 0>(), stage_result.template get<bankSize>()};
	return ret;
}

template<unsigned int NB_CARRY, unsigned int bankSize, unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<QuireDim<N, WES, NB_CARRY>::Size + getNbStages<N, WES, NB_CARRY, bankSize>() - 1, false>
_perform_seg_add_quire(
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
	constexpr unsigned int PROD_EXP_SIZE = PositDim<N, WES>::ProdExpSize;
	constexpr unsigned int PROD_SIGIFICAND_SIZE = PositDim<N, WES>::ProdSignificandSize;
	constexpr unsigned int LOG2_SHIFT_SIZE = Static_Val<getExtShiftSize<N, WES, bankSize>()>::_log2;
	constexpr unsigned int SHIFT_SIZE = getShiftSize<bankSize>();
	constexpr unsigned int EXT_SHIFT_SIZE = getExtShiftSize<N, WES, bankSize>();
	constexpr unsigned int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PROD_SIGIFICAND_SIZE, bankSize>::val;
	constexpr unsigned int PADDING = bankSize*BANKS_FOR_USELESS_BITS - PROD_SIGIFICAND_SIZE;
	constexpr unsigned int EXT_PROD_EXP_SIZE = PROD_EXP_SIZE + BANKS_FOR_USELESS_BITS;

	Wrapper<EXT_PROD_EXP_SIZE, false> padding{{PADDING + 2*PositDim<N, WES>::EMax + 1}};
	auto prod_exp = input.getExp()
						.as_signed()
						.template leftpad<EXT_PROD_EXP_SIZE>()
						.as_unsigned()
						.modularAdd(padding);
	//cerr << to_string(prod_exp) << endl;
	auto input_significand = input.getSignedSignificand();
	auto sign = input.getSignBit();
	auto extended_isSub = Wrapper<PROD_SIGIFICAND_SIZE + 1, false>::generateSequence(isSub);
	auto complementedInputIfIsSub = extended_isSub ^ input_significand;
	//cerr << to_string(complementedInputIfIsSub) << endl;

	auto shift_value = prod_exp.template slice<SHIFT_SIZE - 1, 0>();
	auto ext = complementedInputIfIsSub
							.as_signed()
							.template leftpad<(1<<LOG2_SHIFT_SIZE)>()
							.as_unsigned();
	auto shifted_input = hint::shifter<false, (1<<LOG2_SHIFT_SIZE), SHIFT_SIZE, false>(
			ext,
			shift_value,
			isSub
		);
	auto shifted_input_shrinked = shifted_input.template slice<EXT_SHIFT_SIZE-1, 0>();
	auto stage_select = prod_exp.template slice<PROD_EXP_SIZE-1, SHIFT_SIZE>();

	auto result_is_NaR = input.getIsNaR() | quire.getIsNaR();

	/*
	cerr << "----------------------------------" << endl;
	cerr << to_string(stage_select) << endl;
	cerr << to_string(shifted_input_shrinked) << endl;
	cerr << "=================================" << endl;
	*/
	auto fullQuireWoNaR = _perform_seg_add_quire(
			quire,
			stage_select,
			sign,
			isSub,
			shifted_input_shrinked
		);
	return SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize>{
			result_is_NaR.concatenate(fullQuireWoNaR)
		};
}

template<unsigned int NB_CARRY, unsigned int bankSize, unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int loop_idx>
inline SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize>
_propag_carry_seg_req(
		SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> quire,
		typename enable_if<(loop_idx > 0)>::type* = 0
		)
{
	constexpr unsigned int sssize = PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>();
	constexpr unsigned int essize = getExtShiftSize<N, WES, bankSize>();

	SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> res {
			quire.getIsNaR().concatenate(
				_perform_seg_add_quire(quire,
									   Wrapper<sssize, false>{{0}},
									   Wrapper<1, false>{{0}},
									   Wrapper<1, false>{{0}},
									   Wrapper<essize, false>{{0}}
									)
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
	constexpr unsigned int sssize = PositDim<N, WES>::ProdExpSize - getShiftSize<bankSize>();
	constexpr unsigned int essize = getExtShiftSize<N, WES, bankSize>();

	SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> res {
			quire.getIsNaR().concatenate(
				_perform_seg_add_quire(quire,
									   Wrapper<sssize, false>{{0}},
									   Wrapper<1, false>{{0}},
									   Wrapper<1, false>{{0}},
									   Wrapper<essize, false>{{0}}
									)
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
