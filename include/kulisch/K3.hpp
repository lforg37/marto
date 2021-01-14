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


template<unsigned int WE, unsigned int WF, unsigned int bankSize, unsigned int stage, template<unsigned int, bool> class Wrapper>
Wrapper<bankSize*(stage+1), false> concatAccBanksRec(array<Wrapper<bankSize, false>, getNbStages<WE, WF, bankSize>()> const & acc,
	typename enable_if<(stage == 0)>::type* = 0
	){
	return acc[stage];
}

template<unsigned int WE, unsigned int WF, unsigned int bankSize, unsigned int stage, template<unsigned int, bool> class Wrapper>
Wrapper<bankSize*(stage+1), false> concatAccBanksRec(array<Wrapper<bankSize, false>, getNbStages<WE, WF, bankSize>()> const & acc,
	typename enable_if<(stage >= 1)>::type* = 0
	){
	auto part = concatAccBanksRec<WE, WF, bankSize, stage-1>(acc);
	auto top = acc[stage];
	auto res = top.concatenate(part);
	return res;
}

template<unsigned int WE, unsigned int WF, int unsigned bankSize, template<unsigned int, bool> class Wrapper>
KulischAcc<WE, WF, Wrapper> concatAccBanks(array<Wrapper<bankSize, false>, getNbStages<WE, WF, bankSize>()> const & acc){
	auto res = concatAccBanksRec<WE, WF, bankSize, getNbStages<WE, WF, bankSize>()-1>(acc);
	return res.template slice<IEEEDim<WE, WF>::ACC_SIZE-1, 0>();
}

template<unsigned int WE, unsigned int WF, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
class acc_2CK3
{
		static_assert (hint::is2Pow<bankSize>(), "bankSize should be a power of two");
	private:
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

		array<Wrapper<bankSize, false>, NB_STAGES> banks;
		array<Wrapper<1, false>, NB_STAGES> carries;


		template<unsigned int index>
		inline void setBank(Wrapper<bankSize, false> bank){
			static_assert(index < NB_STAGES, "Trying to set an unexisting bank");
			banks[index] = bank;
		}

		template<unsigned int index>
		inline void setCarry(Wrapper<1, false> carry){
			static_assert (index < NB_STAGES, "Invalid carry index");
			carries[index] = carry;
		}

		template<unsigned int level, unsigned int spread>
		struct IsValidSpreadConfig
		{
			static constexpr bool check = ((level - (spread -1)) <= MAX_STAGE_ID ) and ((static_cast<int>(level) - static_cast<int>(spread) + 1) >= 0);
		};

		template<unsigned int level>
		static constexpr bool canSignExt()
		{
			return (level >= SIGNIFICAND_SPREAD);
		}

		template <unsigned int stageIndex, unsigned int spread>
		inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
			Wrapper<STAGE_ID_WIDTH, false> const & stageSelect,
			Wrapper<1, false> const & inputSign,
			Wrapper<SHIFTED_SIGNIF_WIDTH, false> const &,
			typename enable_if<(spread == 0)>::type* = 0,
			typename enable_if<canSignExt<stageIndex>()>::type* = 0
		) const
		{
			auto cond = stageSelect  < Wrapper<STAGE_ID_WIDTH, false>{{stageIndex+1-SIGNIFICAND_SPREAD}};
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

		template <unsigned int stageIndex, unsigned int spread>
		inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
			Wrapper<STAGE_ID_WIDTH, false> const &,
			Wrapper<1, false> const &,
			Wrapper<SHIFTED_SIGNIF_WIDTH, false> const &,
			typename enable_if<(spread == 0)>::type* = 0,
			typename enable_if<!canSignExt<stageIndex>()>::type* = 0
		) const
		{
			auto ret = Wrapper<bankSize, false>{{0}};
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> add_2CK3_to_add_Rec<cannot_sign_ext, " << stageIndex << ", " << spread << ",...>" << endl;
			cerr << "ret: " << to_string(ret) << endl;
			cerr << "<<< add_2CK3_to_add_Rec<cannot_sign_ext, " << stageIndex << ", " << spread << ",...>" << endl;
		#endif
			return ret;
		}

		template <unsigned int stageIndex, unsigned int spread>
		inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
			Wrapper<STAGE_ID_WIDTH, false> const & stageSelect,
			Wrapper<1, false> const & inputSign,
			Wrapper<SHIFTED_SIGNIF_WIDTH, false> const & shiftedSignificand,
			typename enable_if<(spread >= 1)>::type* = 0,
			typename enable_if<!IsValidSpreadConfig<stageIndex, spread>::check>::type* = 0
			) const
		{
			auto res = add_2CK3_to_add_Rec<stageIndex, spread-1>(stageSelect, inputSign, shiftedSignificand);
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> add_2CK3_to_add_Rec<invalid_config, " << stageIndex << ", " << spread << ",...>" << endl;
			cerr << "res: " << to_string(res) << endl;
			cerr << "<<< add_2CK3_to_add_Rec<invalid_config, " << stageIndex << ", " << spread << ",...>" << endl;
		#endif
			return res;
		}

		template <unsigned int stageIndex, unsigned int spread>
		inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
			Wrapper<STAGE_ID_WIDTH, false> const & stageSelect,
			Wrapper<1, false> const & inputSign,
			Wrapper<SHIFTED_SIGNIF_WIDTH, false> const & shiftedSignificand,
			typename enable_if<(spread >= 1)>::type* = 0,
			typename enable_if<IsValidSpreadConfig<stageIndex, spread>::check>::type* = 0,
			typename enable_if<(spread < SIGNIFICAND_SPREAD)>::type* = 0
			) const
		{
			constexpr unsigned int ident_x = stageIndex - spread + 1;
			constexpr unsigned int up = spread * bankSize - 1;
			constexpr unsigned int down = (spread - 1)*bankSize;

			Wrapper<STAGE_ID_WIDTH, false> activationValue{ident_x};

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

		template <unsigned int stageIndex, unsigned int spread>
		inline Wrapper<bankSize, false> add_2CK3_to_add_Rec(
			Wrapper<STAGE_ID_WIDTH, false> const & stageSelect,
			Wrapper<1, false> const & inputSign,
			Wrapper<SHIFTED_SIGNIF_WIDTH, false> const & shiftedSignificand,
			typename enable_if<(spread >= 1)>::type* = 0,
			typename enable_if<IsValidSpreadConfig<stageIndex, spread>::check>::type* = 0,
			typename enable_if<(spread == SIGNIFICAND_SPREAD)>::type* = 0
			) const
		{
			constexpr unsigned int ident_x = stageIndex - spread + 1;
			constexpr unsigned int up = SHIFTED_SIGNIF_WIDTH - 1;
			constexpr unsigned int down = (spread - 1)*bankSize;
			constexpr unsigned int pad_width = bankSize - (up - down + 1);

			auto ext = Wrapper<pad_width, false>::generateSequence(inputSign);
			auto low = shiftedSignificand.template slice<up, down>();

			Wrapper<STAGE_ID_WIDTH, false> activationValue{ident_x};

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

		template <unsigned int stageIndex>
		inline Wrapper<bankSize, false> add_2CK3_to_add_(
			Wrapper<STAGE_ID_WIDTH, false> const & stageSelect,
			Wrapper<1, false> const & inputSign,
			Wrapper<SHIFTED_SIGNIF_WIDTH, false> const & shiftedSignificand
				) const {
			return add_2CK3_to_add_Rec<stageIndex, SIGNIFICAND_SPREAD>(stageSelect, inputSign, shiftedSignificand);
		}

		template<unsigned int stageIndex>
		inline Wrapper<bankSize+1, false> add_2CK3_acc_stage(
							Wrapper<STAGE_ID_WIDTH, false> const & stageSelect,
							Wrapper<SHIFTED_SIGNIF_WIDTH, false> const & shiftedSignificand,
							Wrapper<1, false> const & sign,
							Wrapper<1, false> const & = {0},
							typename enable_if<(stageIndex>0)>::type* = 0
				) const
		{
			auto bank = getBank<stageIndex>();
			auto accCarry = getCarry<stageIndex-1>();

			auto toAdd = add_2CK3_to_add_<stageIndex>(stageSelect, sign, shiftedSignificand);
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> add_2CK3_acc_stage<" << stageIndex << ">,...>" << endl;
			cerr << "bank: " << to_string(bank) << endl;
			cerr << "accCarry: " << to_string(accCarry) << endl;
			cerr << "toAdd: " << to_string(toAdd) << endl;
			cerr << "<<< add_2CK3_acc_stage<" << stageIndex << ">,...>" << endl;
		#endif
			return bank.addWithCarry(toAdd, accCarry);
		}

		template<unsigned int stageIndex>
		inline Wrapper<bankSize+1, false> add_2CK3_acc_stage(
							Wrapper<STAGE_ID_WIDTH, false> const & stageSelect,
							Wrapper<SHIFTED_SIGNIF_WIDTH, false> const & shiftedSignificand,
							Wrapper<1, false> const & sign,
							Wrapper<1, false> const & carry0 = {0},
							typename enable_if<(stageIndex==0)>::type* = 0
				) const
		{
			auto bank = getBank<stageIndex>();
			auto accCarry = carry0;
			auto toAdd = add_2CK3_to_add_<stageIndex>(stageSelect, sign, shiftedSignificand);
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> add_2CK3_acc_stage<" << stageIndex << ">,...>" << endl;
			cerr << "bank: " << to_string(bank) << endl;
			cerr << "accCarry: " << to_string(accCarry) << endl;
			cerr << "toAdd: " << to_string(toAdd) << endl;
			cerr << "<<< add_2CK3_acc_stage<" << stageIndex << ">,...>" << endl;
		#endif
			return bank.addWithCarry(toAdd, accCarry);
		}


		template<unsigned int stage>
		inline void _perform_addition_all_stages(
					Wrapper<STAGE_ID_WIDTH, false> const & stageSelect,
					Wrapper<SHIFTED_SIGNIF_WIDTH, false> const & shiftedSignificand,
					Wrapper<1, false> const & sign,
					acc_2CK3& resAcc,
					typename enable_if<(stage == 0)>::type* = 0
				) const
		{
			auto stageResult = add_2CK3_acc_stage<stage>(
						stageSelect,
						shiftedSignificand,
						sign
			);
			resAcc.template setCarry<stage>(stageResult.template get<bankSize>());
			resAcc.template setBank<stage>(stageResult.template slice<bankSize-1,0>());
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> _perform_addition_all_stages<" << stage << ", ...>" << endl;
			cerr << "stageResult: " << to_string(stageResult) << endl;
			cerr << "<<< _perform_addition_all_stages<" << stage << ", ...>" << endl;
		#endif
		}

		/**
		 * Recursive unrolling of addition on all accumulator banks
		 * @param stageSelect[in] range in which the product should be added
		 * @param shiftedSignificand[in] significand shifted by a value of the range (0, banksize(
		 * @param sign[in] sign of the operation (1 for subtraction)
		 * @param resAcc[out] accumulator in which to store the result
		 */
		template<unsigned int stage>
		inline void _perform_addition_all_stages(
					Wrapper<STAGE_ID_WIDTH, false> const & stageSelect,
					Wrapper<SHIFTED_SIGNIF_WIDTH, false> const & shiftedSignificand,
					Wrapper<1, false> const & sign,
					acc_2CK3& resAcc,
					typename enable_if<(stage > 0)>::type* = 0
				) const
		{
			auto stageResult = add_2CK3_acc_stage<stage>(
						stageSelect,
						shiftedSignificand,
						sign
			);
			resAcc.template setCarry<stage>(stageResult.template get<bankSize>());
			resAcc.template setBank<stage>(stageResult.template slice<bankSize-1,0>());
		#ifdef K3_ACC_DEBUG_FULL
			cerr << ">>> _perform_addition_all_stages<" << stage << ", ...>" << endl;
			cerr << "stageResult: " << to_string(stageResult) << endl;
			cerr << "<<< _perform_addition_all_stages<" << stage << ", ...>" << endl;
		#endif
			_perform_addition_all_stages<stage-1>(stageSelect, shiftedSignificand, sign, resAcc);
		}

		template<unsigned int nb_repeat>
		inline acc_2CK3 _propagate_carries_2CK3_Req(
				typename enable_if<(nb_repeat > 0)>::type* = 0
		) const
		{
				acc_2CK3 ret;
				_perform_addition_all_stages<NB_STAGES - 1>({{0}}, {{0}}, {{0}}, ret);
				return ret._propagate_carries_2CK3_Req<nb_repeat - 1>();
		}

		template<unsigned int nb_repeat>
		inline acc_2CK3 _propagate_carries_2CK3_Req(
				typename enable_if<(nb_repeat == 0)>::type* = 0
		) const
		{
				acc_2CK3 ret;
				_perform_addition_all_stages<NB_STAGES - 1>({{0}}, {{0}}, {{0}}, ret);
				return ret;
		}

	public:
		acc_2CK3(
				KulischAcc<WE, WF, Wrapper> const & acc = {{0}},
				Wrapper<NB_STAGES, false> const & carries_in= {{0}}
				)
		{
			ArraySplitter<ACC_SIZE, bankSize>::distribute(
						static_cast<Wrapper<ACC_SIZE, false> const &>(acc), banks);
			ArraySplitter<NB_STAGES, 1>::distribute(carries_in, carries);
#ifdef K3_ACC_DEBUG_FULL
			cerr << "CTOR of acc_2CK3<WE=" << WE << ", WF=" << WF << ", bankSize=" << bankSize << ">" << endl;
			cerr << "ACC_SIZE=" << ACC_SIZE << endl;
			cerr << "NB_STAGES=" << NB_STAGES << endl;
#endif
		}

		inline KulischAcc<WE, WF, Wrapper> getAcc(){
			return concatAccBanks<WE, WF, bankSize>(banks);
		}

		template<unsigned int index>
		inline Wrapper<bankSize, false> getBank() const{
			static_assert(index < NB_STAGES, "Trying to access an unexisting bank");
			return banks[index];
		}

		template<unsigned int index>
		inline Wrapper<1, false> getCarry() const{
			static_assert (index < NB_STAGES, "Invalid carry index");
			return carries[index];
		}

		inline Wrapper<1, false> isNeg() const {
			return banks[NB_STAGES-1].template get<bankSize-1>();
		}

		/**
		 *	Adding a value to a 2's complement segmented Kulisch accumulator
		 *	@param acc [in] the accumulator to which the product should be added
		 *	@param prod [in] the product to add
		 */
		inline acc_2CK3 add_2CK3(FPProd<WE, WF, Wrapper> const & prod) const
		{
			auto prodExp = prod.getExp();
			auto prodSign = prod.getSignBit();

			auto inputSignificand = prod.getSignificand();
			auto inputSignificandComplemented = Wrapper<IEEEDim<WE, WF>::WFF_Prod, false>::mux(
						prodSign,
						(inputSignificand.invert().modularAdd(Wrapper<IEEEDim<WE, WF>::WFF_Prod, false>{1})),
						 inputSignificand
				);

			auto shiftValue = prodExp.template slice<BANK_WIDTH-1,0>();
			auto extension = Wrapper<SHIFTED_SIGNIF_WIDTH - IEEEDim<WE, WF>::WFF_Prod, false>::generateSequence(prodSign);
			auto ext = extension.concatenate(inputSignificandComplemented);

			auto shiftedInput = shifter<false>(ext,shiftValue,{0});

			auto stageSelect = prodExp.template slice<WE, BANK_WIDTH>();
			acc_2CK3<WE, WF, bankSize, Wrapper> fullAcc{};

			_perform_addition_all_stages<NB_STAGES-1>(stageSelect, shiftedInput, prodSign, fullAcc);
		#ifdef K3_ACC_DEBUG
			cerr << "====== K3 accumulation ====" << endl;
			cerr << "prodExp: " << to_string(prodExp) << endl;
			cerr << "prodSign: " << to_string(prodSign) << endl;
			cerr << "inputSignificand: " << to_string(inputSignificand) << endl;
			cerr << "inputSignificandComplemented: " << to_string(inputSignificandComplemented) << endl;
			cerr << "extension: " << to_string(extension) << endl;
			cerr << "shiftValue: " << to_string(shiftValue) << endl;
			cerr << "ext: " << to_string(ext) << endl;
			cerr << "shiftedInput: " << to_string(shiftedInput) << endl;
			cerr << "stageSelect: " << to_string(stageSelect) << endl;
			cerr << "Resulting acc :" << to_string(fullAcc.getAcc().downcast()) << endl;
			size_t i = 0;
			for (auto const & bank : fullAcc.banks) {
				cerr << "banks[" << i++ << "]: " << to_string(bank) << endl;
			}
			cerr << "===========================" << endl;
		#endif
			return fullAcc;
		}

		inline KulischAcc<WE, WF, Wrapper> propagate_carries_2CK3() const
		{
			auto retK3 = _propagate_carries_2CK3_Req<NB_STAGES>();
			auto ret = retK3.getAcc();
		#ifdef K3_ACC_DEBUG
			cerr << "ret: " << to_string(static_cast<Wrapper<FPDim<WE, WF>::ACC_SIZE, false> const &>(ret)) << endl;
		#endif
			return ret;
		}
};
#endif
