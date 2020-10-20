#ifndef IEEETYPE_HPP
#define IEEETYPE_HPP

#include <cstdint>

enum struct IEEERoundingMode : uint8_t {
	RoundTowardZero = 0,
	RoundUp = 7,
	RoundDown = 5,
	RoundNearestTieEven = 4,
	RoundNearestTieAway = 6
};

template<unsigned int N>
struct StandardFPDim;

template<unsigned int _WE, unsigned int _WF>
struct FPDim {
		static constexpr unsigned int WE = _WE;
		static constexpr unsigned int WF = _WF;
		static constexpr unsigned int BIAS = (1<<(WE-1))-1;
		static constexpr unsigned int ACC_SIZE = (1<<(WE+1)) -1 + 2*WF+2;
		static constexpr unsigned int ACC_MID = (1<<(WE)) -1 + 2*WF+2;
		static constexpr unsigned int FP_SPREAD = (1<<(WE-1));
		static constexpr unsigned int PROD_FP_SPREAD = (1<<(WE));
		static constexpr unsigned int SUBNORMAL_LIMIT = (1<<(WE))+BIAS;
};

template<>
struct StandardFPDim<16> : public FPDim<5, 10> {};

template <>
struct StandardFPDim<32> : public FPDim<8, 23> {};

template<>
struct StandardFPDim<64> : public FPDim<11, 52> {};

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
	class IEEENumber : public Wrapper<WE + WF + 1, false>
	{
		private:
			typedef Wrapper<WE+WF+1, false> basetype;
			template<unsigned int W>
			using us_wrapper = Wrapper<W, false>;
		public:
			static constexpr unsigned int _WE = WE;
			static constexpr unsigned int _WF = WF;
			using rounding_type_t = us_wrapper<3>;

			IEEENumber(Wrapper<WE + WF + 1, false> const & val):Wrapper<WE+WF+1, false>{val}{}

			inline us_wrapper<1> getSign() const{
				return basetype::template get<WE+WF>();
			}

			inline us_wrapper<WE> getExponent() const {
				return basetype::template slice<WF + WE - 1, WF>();
			}

			inline us_wrapper<WF> getFractionnalPart() const {
				return basetype::template slice<WF - 1, 0>();
			}

			inline us_wrapper<WF+WE> getExpFrac() const {
				return basetype::template slice<WF+WE - 1, 0>();
			}

			inline us_wrapper<1> getLeadBitVal() const {
				return getExponent().or_reduction();
			}

			inline us_wrapper<1> isInfinity() const {
				return getExponent().and_reduction().bitwise_and(getFractionnalPart().nor_reduction());
			}

			inline us_wrapper<1> isNaN() const {
				return (getExponent().and_reduction()) & (getFractionnalPart().or_reduction());
			}
	};
#endif // IEEETYPE_HPP
