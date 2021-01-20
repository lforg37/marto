#ifndef K3_HPP
#define K3_HPP
#include "kulisch_dim.hpp"
using namespace std;
#include "primitives/lzoc_shifter.hpp"

using hint::Static_Val;
using hint::Static_Ceil_Div;
using hint::shifter;

#ifdef K3_ACC_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif

template<unsigned int WE, unsigned int WF, int bankSize>
static constexpr unsigned int getNbStages(){
	return Static_Ceil_Div<IEEEDim<WE, WF>::ACC_SIZE, bankSize>::val;
}

template<int N, int bankSize>
static constexpr int getSegmentedAccSize(){
	return getNbStages<N, bankSize>() * bankSize;
}

template<unsigned int WE, unsigned int WF, unsigned int bankSize>
static constexpr unsigned int getMantSpread(){
	return Static_Ceil_Div<IEEEDim<WE, WF>::WFF_Prod, bankSize>::val+1;
}

template<unsigned int WE, unsigned int WF, unsigned int bankSize>
struct AccK3Dim
{
		using _dim = IEEEDim<WE, WF>;
		/**
		 * @brief NB_STAGES required number of banks inside the accumulator
		 */

		static constexpr unsigned int NB_STAGES = getNbStages<WE, WF, bankSize>();
		static constexpr unsigned int ACC_SIZE = _dim::ACC_SIZE;

		/**
		 * @brief BANK_WIDTH size of the index inside the bank
		 */
		static constexpr unsigned int BANK_WIDTH = Static_Val<bankSize - 1>::_storage;

		/**
		 * @brief SHIFTED_SIGNIF_WIDTH Width of the shifted significand.
		 */
		static constexpr unsigned int SHIFTED_SIGNIF_WIDTH = _dim::WFF_Prod + bankSize - 1;

		/**
		 * @brief STAGE_ID_WIDTH size of the exponent part which identifies the first stage on which the add should start.
		 */
		static constexpr unsigned int STAGE_ID_WIDTH = _dim::WE_Prod - BANK_WIDTH;

		/**
		 * @brief SIGNIFICAND_SPREAD Maximum number of banks over which the shifted significand can spread
		 */
		static constexpr unsigned int SIGNIFICAND_SPREAD = getMantSpread<WE, WF, bankSize>();

		static constexpr unsigned int MAX_STAGE_ID = (1 << STAGE_ID_WIDTH) - 1;

		static constexpr unsigned int TOTAL_WIDTH = (bankSize + 1) * NB_STAGES;
};

template<unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
class acc_2CK3 : public Wrapper<AccK3Dim<WE, WF, bankSize>::TOTAL_WIDTH, false>
{
		static_assert (hint::is2Pow<bankSize>(), "bankSize should be a power of two");
	private:

		using _k3dim = AccK3Dim<WE, WF, bankSize>;

		using _storage = Wrapper<_k3dim::TOTAL_WIDTH, false>;

	public:
		acc_2CK3(
				Wrapper<_k3dim::TOTAL_WIDTH, false> const & acc = {{0}}
				):_storage{acc}
		{}

		template<unsigned int index>
		inline Wrapper<bankSize, false> getBank() const{
			static_assert(index < _k3dim::NB_STAGES, "Trying to access an unexisting bank");
			constexpr unsigned int low_idx = index*bankSize;
			constexpr unsigned int up_idx = low_idx + bankSize - 1;
			return _storage::template slice<up_idx, low_idx>();
		}

		template<unsigned int index>
		inline Wrapper<1, false> getCarry() const{
			static_assert (index < _k3dim::NB_STAGES, "Invalid carry index");
			return _storage::template get<index + _k3dim::NB_STAGES * bankSize>();
		}

		inline Wrapper<1, false> isNeg() const {
			return _storage::template get<_k3dim::NB_STAGES * bankSize - 1>();
		}
};

template <unsigned int WE, unsigned int WF, unsigned int bankSize>
class AccK3Adder {
	private:
		using _k3dim = AccK3Dim<WE, WF, bankSize>;

		template <template<unsigned int, bool> class Wrapper>
		using prod_t = FPProd<WE, WF, Wrapper>;

		template <template<unsigned int, bool> class Wrapper>
		using k3acc_t = acc_2CK3<WE, WF, bankSize, Wrapper>;

		template <template<unsigned int, bool> class Wrapper>
		using kulisch_t = KulischAcc<WE, WF, Wrapper>;

		static constexpr unsigned int STAGE_ID_WIDTH = WE-_k3dim::BANK_WIDTH+1;

		template <template<unsigned int, bool> class Wrapper>
		using stage_id_t = Wrapper<STAGE_ID_WIDTH, false>;

		template <template<unsigned int, bool> class Wrapper>
		using shifted_signif_t = Wrapper<_k3dim::SHIFTED_SIGNIF_WIDTH, false>;

		template<unsigned int loop_idx, template<unsigned int, bool> class Wrapper>
		struct _partial_acc_store {
				Wrapper<loop_idx + 1, false> carries;
				Wrapper<(loop_idx + 1)*bankSize, false> banks;
		};

		template<unsigned int level, unsigned int spread>
		struct IsValidSpreadConfig
		{
			static constexpr bool check = ((level - (spread -1)) <= _k3dim::MAX_STAGE_ID ) and ((static_cast<int>(level) - static_cast<int>(spread) + 1) >= 0);
		};

		template<unsigned int level>
		struct CanSignExt
		{
			static constexpr bool check = (level >= _k3dim::SIGNIFICAND_SPREAD);
		};

		template <unsigned int stageIndex, unsigned int spread, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
			stage_id_t<Wrapper> const & stageSelect,
			Wrapper<1, false> const & inputSign,
			shifted_signif_t<Wrapper> const &,
			typename enable_if<(spread == 0)>::type* = 0,
			typename enable_if<CanSignExt<stageIndex>::check>::type* = 0
		)
		{
			auto cond = stageSelect  < Wrapper<_k3dim::STAGE_ID_WIDTH, false>{{stageIndex+1-_k3dim::SIGNIFICAND_SPREAD}};
			auto ret = Wrapper<bankSize, false>::generateSequence(cond & inputSign);
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> add_2CK3_to_add_Rec<can_sign_ext, " << stageIndex << ", " << spread << ",...>" << endl;
			cerr << "stageSelect: " << to_string(stageSelect) << endl;
			cerr << "cond: " << to_string(cond) << endl;
			cerr << "ret: " << to_string(ret) << endl;
			cerr << "<<< add_2CK3_to_add_Rec<can_sign_ext, " << stageIndex << ", " << spread << ",...>" << endl;
		#endif
			return ret;
		}

		template <unsigned int stageIndex, unsigned int spread, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
			stage_id_t<Wrapper> const &,
			Wrapper<1, false> const &,
			shifted_signif_t<Wrapper> const &,
			typename enable_if<(spread == 0)>::type* = 0,
			typename enable_if<!CanSignExt<stageIndex>::check>::type* = 0
		)
		{
			auto ret = Wrapper<bankSize, false>{{0}};
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> add_2CK3_to_add_Rec<cannot_sign_ext, " << stageIndex << ", " << spread << ",...>" << endl;
			cerr << "ret: " << to_string(ret) << endl;
			cerr << "<<< add_2CK3_to_add_Rec<cannot_sign_ext, " << stageIndex << ", " << spread << ",...>" << endl;
		#endif
			return ret;
		}

		template <unsigned int stageIndex, unsigned int spread, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
			stage_id_t<Wrapper> const & stageSelect,
			Wrapper<1, false> const & inputSign,
			shifted_signif_t<Wrapper> const & shiftedSignificand,
			typename enable_if<(spread >= 1)>::type* = 0,
			typename enable_if<!IsValidSpreadConfig<stageIndex, spread>::check>::type* = 0
			)
		{
			auto res = add_2CK3_to_add_Rec<stageIndex, spread-1>(stageSelect, inputSign, shiftedSignificand);
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> add_2CK3_to_add_Rec<invalid_config, " << stageIndex << ", " << spread << ",...>" << endl;
			cerr << "res: " << to_string(res) << endl;
			cerr << "<<< add_2CK3_to_add_Rec<invalid_config, " << stageIndex << ", " << spread << ",...>" << endl;
		#endif
			return res;
		}

		template <unsigned int stageIndex, unsigned int spread, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
			stage_id_t<Wrapper> const & stageSelect,
			Wrapper<1, false> const & inputSign,
			shifted_signif_t<Wrapper> const & shiftedSignificand,
			typename enable_if<(spread >= 1)>::type* = 0,
			typename enable_if<IsValidSpreadConfig<stageIndex, spread>::check>::type* = 0,
			typename enable_if<(spread < _k3dim::SIGNIFICAND_SPREAD)>::type* = 0
			)
		{
			constexpr unsigned int ident_x = stageIndex - spread + 1;
			constexpr unsigned int up = spread * bankSize - 1;
			constexpr unsigned int down = (spread - 1)*bankSize;

			stage_id_t<Wrapper> activationValue{ident_x};

			auto cond = (stageSelect == activationValue);

			auto res = Wrapper<bankSize, false>::mux(
						cond,
						shiftedSignificand.template slice<up, down>(),
						add_2CK3_to_add_Rec<stageIndex, spread-1>(stageSelect, inputSign, shiftedSignificand)
				);
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> add_2CK3_to_add_Rec<valid_config, " << stageIndex << ", " << spread << ",...>" << endl;
			cerr << "stageSelect: " << to_string(stageSelect) << endl;
			cerr << "cond: " << to_string(cond) << endl;
			cerr << "res: " << to_string(res) << endl;
			cerr << "<<< add_2CK3_to_add_Rec<valid_config, " << stageIndex << ", " << spread << ",...>" << endl;
		#endif
			return res;
		}

		template <unsigned int stageIndex, unsigned int spread, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
			stage_id_t<Wrapper> const & stageSelect,
			Wrapper<1, false> const & inputSign,
			shifted_signif_t<Wrapper> const & shiftedSignificand,
			typename enable_if<(spread >= 1)>::type* = 0,
			typename enable_if<IsValidSpreadConfig<stageIndex, spread>::check>::type* = 0,
			typename enable_if<(spread == _k3dim::SIGNIFICAND_SPREAD)>::type* = 0
			)
		{
			constexpr unsigned int ident_x = stageIndex - spread + 1;
			constexpr unsigned int up = _k3dim::SHIFTED_SIGNIF_WIDTH - 1;
			constexpr unsigned int down = (spread - 1)*bankSize;
			constexpr unsigned int pad_width = bankSize - (up - down + 1);

			auto ext = Wrapper<pad_width, false>::generateSequence(inputSign);
			auto low = shiftedSignificand.template slice<up, down>();

			stage_id_t<Wrapper> activationValue{ident_x};

			auto cond = (stageSelect == activationValue);

			auto res = Wrapper<bankSize, false>::mux(
						cond,
						ext.concatenate(low),
						add_2CK3_to_add_Rec<stageIndex, spread-1>(stageSelect, inputSign, shiftedSignificand)
				);
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> add_2CK3_to_add_Rec<valid_config, " << stageIndex << ", " << spread << ",...>" << endl;
			cerr << "stageSelect: " << to_string(stageSelect) << endl;
			cerr << "cond: " << to_string(cond) << endl;
			cerr << "res: " << to_string(res) << endl;
			cerr << "<<< add_2CK3_to_add_Rec<valid_config, " << stageIndex << ", " << spread << ",...>" << endl;
		#endif
			return res;
		}

		template <unsigned int stageIndex, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<bankSize, false> _get_to_add(
			stage_id_t<Wrapper> const & stageSelect,
			Wrapper<1, false> const & inputSign,
			shifted_signif_t<Wrapper> const & shiftedSignificand
				)  {
			return add_2CK3_to_add_Rec<stageIndex, _k3dim::SIGNIFICAND_SPREAD>(stageSelect, inputSign, shiftedSignificand);
		}



		template<unsigned int loop_idx, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<bankSize+1, false> _add_sub_acc_stage(
				k3acc_t<Wrapper> const & acc,
				stage_id_t<Wrapper> const & stageSelect,
				Wrapper<1, false> const & effectSub,
				shifted_signif_t<Wrapper> const & shiftedInput,
				typename enable_if<(loop_idx > 0)>::type* = 0
			)
		{
			auto accBank = acc.template getBank<loop_idx>();
			auto carry = acc.template getCarry<loop_idx - 1>();

			auto toAdd = _get_to_add<loop_idx>(stageSelect, effectSub, shiftedInput);

			auto sum = accBank.addWithCarry(toAdd, carry);
			return sum;
		}


		template<unsigned int loop_idx, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<bankSize+1, false> _add_sub_acc_stage(
				k3acc_t<Wrapper> const & acc,
				stage_id_t<Wrapper> const & stageSelect,
				Wrapper<1, false> const & effectSub,
				shifted_signif_t<Wrapper> const & shiftedInput,
				typename enable_if<(loop_idx == 0)>::type* = 0
			)
		{
			auto accBank = acc.template getBank<loop_idx>();

			auto toAdd = _get_to_add<loop_idx>(stageSelect, effectSub, shiftedInput);

			auto sum = accBank + toAdd;
			return sum;
		}

		template<unsigned int loop_idx, template<unsigned int, bool> class Wrapper>
		static inline _partial_acc_store<loop_idx, Wrapper> _req_partial_add(
					k3acc_t<Wrapper> const & acc,
					stage_id_t<Wrapper> stageSelect,
					Wrapper<1, false> effectSub,
					shifted_signif_t<Wrapper> shiftedInput,
					typename enable_if<(loop_idx > 0)>::type* = 0
				)
		{
			auto stage_result = _add_sub_acc_stage<loop_idx>(
						acc,
						stageSelect,
						effectSub,
						shiftedInput
				);
			auto lower_stages = _req_partial_add<loop_idx - 1>(acc, stageSelect, effectSub, shiftedInput);
			auto cur_stage_carry = stage_result.template get<bankSize>();
			auto cur_stage_bank = stage_result.template slice<bankSize - 1, 0>();
			_partial_acc_store<loop_idx, Wrapper> ret {
					cur_stage_carry.concatenate(lower_stages.carries),
					cur_stage_bank.concatenate(lower_stages.banks)
			};
			return ret;
		}

		template<unsigned int loop_idx, template<unsigned int, bool> class Wrapper>
		static inline _partial_acc_store<loop_idx, Wrapper> _req_partial_add(
					k3acc_t<Wrapper> const & acc,
					stage_id_t<Wrapper> stageSelect,
					Wrapper<1, false> effectSub,
					shifted_signif_t<Wrapper> shiftedInput,
					typename enable_if<(loop_idx== 0)>::type* = 0
				)
		{
			 auto stage_result = _add_sub_acc_stage<0>(
						 acc,
						 stageSelect,
						 effectSub,
						 shiftedInput
				 );
			_partial_acc_store<loop_idx, Wrapper> ret{stage_result.template get<bankSize>(), stage_result.template slice<bankSize - 1, 0>()};
			return ret;
		}

		template<template<unsigned int, bool> class Wrapper>
		static inline k3acc_t<Wrapper>_perform_seg_add(
				k3acc_t<Wrapper> const & acc,
				stage_id_t<Wrapper> const & stageSelect,
				Wrapper<1, false> const & effectSub,
				shifted_signif_t<Wrapper> const & shiftedInput)
		{
			constexpr unsigned int start_idx = _k3dim::NB_STAGES - 1;
			auto val = _req_partial_add<start_idx>(
					acc,
					stageSelect,
					effectSub,
					shiftedInput
				);
			return {{val.carries.concatenate(val.banks)}};
		}


	public:
		template<template<unsigned int, bool> class Wrapper>
		static inline k3acc_t<Wrapper> add_2CK3(
				k3acc_t<Wrapper> const & in,
				prod_t<Wrapper> const & prod,
				Wrapper<1, false> const & isSub)
		{
			auto prodExp = prod.getExp();
			auto prodSign = prod.getSignBit();
			auto effectSub = prodSign ^ isSub;

			auto inputSignificand = prod.getSignificand();
			auto inputSignificandComplemented = Wrapper<IEEEDim<WE, WF>::WFF_Prod, false>::mux(
						effectSub,
						(inputSignificand.invert().modularAdd(Wrapper<IEEEDim<WE, WF>::WFF_Prod, false>{1})),
						 inputSignificand
				);

			auto shiftValue = prodExp.template slice<_k3dim::BANK_WIDTH-1,0>();
			auto extension = Wrapper<_k3dim::SHIFTED_SIGNIF_WIDTH - IEEEDim<WE, WF>::WFF_Prod, false>::generateSequence(prodSign);
			auto ext = extension.concatenate(inputSignificandComplemented);

			auto shiftedInput = shifter<false>(ext,shiftValue,{0});

			auto stageSelect = prodExp.template slice<WE, _k3dim::BANK_WIDTH>();

			auto res = _perform_seg_add(in, stageSelect, effectSub, shiftedInput);
			return res;
		}
};

template<unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline acc_2CK3<WE, WF, bankSize, Wrapper> add_2ck3(
		acc_2CK3<WE, WF, bankSize, Wrapper> const & acc,
		FPProd<WE, WF, Wrapper> const & prod,
		Wrapper<1, false> const & issub = {{0}}
	)
{
	return AccK3Adder<WE, WF, bankSize>::add_2CK3(acc, prod, issub);
}

template <unsigned int WE, unsigned int WF, unsigned int bankSize>
class K3ToKulisch {
	private:
		using _k3dim = AccK3Dim<WE, WF, bankSize>;

		template <template<unsigned int, bool> class Wrapper>
		using prod_t = FPProd<WE, WF, Wrapper>;

		template <template<unsigned int, bool> class Wrapper>
		using k3acc_t = acc_2CK3<WE, WF, bankSize, Wrapper>;

		template <template<unsigned int, bool> class Wrapper>
		using kulisch_t = KulischAcc<WE, WF, Wrapper>;

		static constexpr unsigned int NB_STAGES = _k3dim::NB_STAGES;

		using adder = AccK3Adder<WE, WF, bankSize>;


		template<unsigned int loop_idx, template<unsigned int, bool> class Wrapper>
		static inline kulisch_t<Wrapper> _propagate_rec(
			k3acc_t<Wrapper> const & acc,
			typename enable_if<(loop_idx > 0)>::type* = 0
		)
		{
			auto next_acc = adder::add_2CK3(acc, {{0}}, {{0}});
			return _propagate_rec<loop_idx - 1>(next_acc);
		}

		template<unsigned int loop_idx, template<unsigned int, bool> class Wrapper>
		static inline kulisch_t<Wrapper> _propagate_rec(
			k3acc_t<Wrapper> const & acc,
			typename enable_if<(loop_idx == 0)>::type* = 0
		)
		{
			auto next_acc = adder::add_2CK3(acc, {{0}}, {{0}});
			return next_acc.template slice<_k3dim::ACC_SIZE-1, 0>();
		};

	public:
		template<template<unsigned int, bool> class Wrapper>
		static inline kulisch_t<Wrapper> propagate_carries(k3acc_t<Wrapper> const & in)
		{
			return _propagate_rec<NB_STAGES-1>(in);
		}
};

template<unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
inline KulischAcc<WE, WF, Wrapper> propagate_carries(acc_2CK3<WE, WF, bankSize, Wrapper> const & acc)
{
	return K3ToKulisch<WE, WF, bankSize>::propagate_carries(acc);
}
#endif
