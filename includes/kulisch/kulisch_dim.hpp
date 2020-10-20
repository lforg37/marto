#ifndef FP_DIM_TPP
#define FP_DIM_TPP

#include <cstdint>
#include <type_traits>

#include "tools/static_math.hpp"
#include "helpers/splitting.hpp"
#include "ieeefloats/ieeetype.hpp"

using namespace std;

using hint::Static_Val;
using hint::Static_Ceil_Div;
using hint::ArraySplitter;

template<unsigned int WE_, unsigned int WF_>
struct IEEEProdDim
{
	static constexpr unsigned int WProdExp = WE_ + 1;
	static constexpr unsigned int WProdFrac = 2*WF_ + 2;
	static constexpr unsigned int WProd = 1 + WProdExp + WProdFrac;
};

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
class FPProd : public Wrapper<IEEEProdDim<WE, WF>::WProd, false>
{
	private:
		constexpr static unsigned int _width = IEEEProdDim<WE, WF>::WProd;
		using storage_type = Wrapper<_width, false>;
		using _dim = IEEEProdDim<WE, WF>;
	public:
		FPProd(
				Wrapper<1, false> mult_s,
				Wrapper<WE+1, false> mult_e,
				Wrapper<2*WF+2, false> mult_m
		)
		{
			auto signed_exp = mult_s.concatenate(mult_e);
			storage_type{signed_exp.concatenate(mult_m)};
		}

		FPProd(storage_type const & val):storage_type{val}{}

		inline Wrapper<_dim::WProdFrac, false> getSignificand()
		{
			return storage_type::template slice<_dim::WProdFrac - 1, 0>();
		}


		inline Wrapper<1, false> getSignBit()
		{
			return storage_type::template get<_dim::WProd-1>();
		}

		inline Wrapper<_dim::WProdExp, false> getExp()
		{
			return storage_type::template slice<_dim::WProd-2, _dim::WProdFrac>();
		}
};

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
using KulischAcc = Wrapper<FPDim<WE, WF>::ACC_SIZE, false>;

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
using SignedKulischAcc = Wrapper<FPDim<WE, WF>::ACC_SIZE+1, false>;


template<unsigned int WE, unsigned int WF, int bankSize>
static constexpr unsigned int getNbStages(){
	return Static_Ceil_Div<FPDim<WE, WF>::ACC_SIZE, bankSize>::val;
}

template<int N, int bankSize>
static constexpr int getSegmentedAccSize(){
	return getNbStages<N, bankSize>() * bankSize;
}

template<unsigned int WE, unsigned int WF, unsigned int bankSize>
static constexpr unsigned int getMantSpread(){
	return Static_Ceil_Div<IEEEProdDim<WE, WF>::WProdFrac,bankSize>::val+1;
}


template<unsigned int WE, unsigned int WF, unsigned int bankSize, unsigned int stage, template<unsigned int, bool> class Wrapper>
Wrapper<bankSize*(stage+1), false> concatAccBanksRec(Wrapper<bankSize, false> const acc[getNbStages<WE, WF, bankSize>()],
	typename enable_if<(stage == 0)>::type* = 0
	){
	return acc[stage];
}

template<unsigned int WE, unsigned int WF, int bankSize, int stage, template<unsigned int, bool> class Wrapper>
Wrapper<bankSize*(stage+1), false> concatAccBanksRec(Wrapper<bankSize, false> const acc[getNbStages<WE, WF, bankSize>()],
	typename enable_if<(stage >= 1)>::type* = 0
	){
	auto part = concatAccBanksRec<WE, WF, bankSize, stage-1>(acc);
	auto top = acc[stage];
	auto res = top.concat(part);
	return res;
}

template<unsigned int WE, unsigned int WF, int bankSize, template<unsigned int, bool> class Wrapper>
KulischAcc<WE, WF, Wrapper> concatAccBanks(Wrapper<bankSize, false> acc[getNbStages<WE, WF, bankSize>()]){
	auto res = concatAccBanksRec<WE, WF, bankSize, getNbStages<WE, WF, bankSize>()-1>(acc);
	return res.template slice<FPDim<WE, WF>::ACC_SIZE-1, 0>();
}



template<unsigned int WE, unsigned int WF, int bankSize, template<unsigned int, bool> class Wrapper>
class acc_2CK3
{
	private:
		static constexpr unsigned int NB_STAGES = getNbStages<WE, WF, bankSize>();

		Wrapper<bankSize, false> banks[NB_STAGES];
		Wrapper<1, false> carries[NB_STAGES];



	public:
		acc_2CK3(
				KulischAcc<WE, WF, Wrapper> acc = {{0}},
				Wrapper<NB_STAGES, false> carries_in = {{0}}
				)
		{
			ArraySplitter<FPDim<WE, WF>::ACC_SIZE, bankSize>::distribute(acc, banks);
			ArraySplitter<NB_STAGES, 1>::distribute(carries_in, carries);
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
		inline void setBank(Wrapper<bankSize, false> bank){
			static_assert(index < NB_STAGES, "Trying to set an unexisting bank");
			banks[index] = bank;
		}

		template<unsigned int index>
		inline Wrapper<1, false> getCarry() const{
			static_assert (index < NB_STAGES, "Invalid carry index");
			return carries[index];
		}

		template<unsigned int index>
		inline void setCarry(Wrapper<1, false> carry){
			static_assert (index < NB_STAGES, "Invalid carry index");
			carries[index] = carry;
		}

		inline Wrapper<1, false> isNeg() const {
			return banks[NB_STAGES-1].template get<bankSize-1>();
		}
};
#endif
