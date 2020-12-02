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
struct StandardIEEEDim;

template<unsigned int _WE, unsigned int _WF>
struct IEEEDim {
		/**
		 * @brief Format width
		 */
		static constexpr unsigned int W = _WE + _WF + 1;

		/**
		 * @brief WE Exponent field width
		 */
		static constexpr unsigned int WE = _WE;

		/**
		 * @brief WF Fractional width
		 */
		static constexpr unsigned int WF = _WF;

		/**
		 * @brief BIAS IEEE-754 bias value
		 */
		static constexpr unsigned int BIAS = (1<<(WE-1))-1;

		/**
		 * @brief WFF Full significand width
		 */
		static constexpr unsigned int WFF = WF +1;

		/**
		 * @brief MAX_NORMAL_BIASED_EXP Maximum valid exponent in biased form
		 */
		static constexpr unsigned int MAX_NORMAL_BIASED_EXP = (1 << WE) - 2;

		/**
		 * @brief MIN_NORMAL_BIASED_EXP Minimal normal number exponent in biased form
		 */
		static constexpr unsigned int MIN_NORMAL_BIASED_EXP = 1;

		/**
		 * @brief DIFF_EXP_NB Number of different valid exponent values
		 */
		static constexpr unsigned int DIFF_EXP_NB = MAX_NORMAL_BIASED_EXP - MIN_NORMAL_BIASED_EXP + 1;

		/**
		 * @brief MAX_NORMAL_UNBIASED_EXP Maximum valid exponent of the format
		 */
		static constexpr int MAX_NORMAL_UNBIASED_EXP = static_cast<int>(MAX_NORMAL_BIASED_EXP) - static_cast<int>(BIAS);

		/**
		 * @brief MIN_NORMAL_UNBIASED_EXP Minimum valid exponent of the format
		 */
		static constexpr int MIN_NORMAL_UNBIASED_EXP = static_cast<int>(MIN_NORMAL_BIASED_EXP) - static_cast<int>(BIAS);

		/**
		 * @brief MIN_LOWBIT_WEIGHT Minimal weight of the last fraction bit
		 */
		static constexpr int MIN_LOWBIT_WEIGHT = MIN_NORMAL_UNBIASED_EXP - static_cast<int>(WF);

		/**
		 * @brief MEX_LOWBIT_WEIGHT Minimal weight of the last fraction bit
		 */
		static constexpr int MAX_LOWBIT_WEIGHT = MAX_NORMAL_UNBIASED_EXP - static_cast<int>(WF);

		/**
		 * @brief WE_Prod two operand product exponent field width
		 */
		static constexpr unsigned int WE_Prod = WE+1;

		/**
		 * @brief WFF_Prod two op exact product full significand width
		 */
		static constexpr unsigned int WFF_Prod = 2*WFF;

		/**
		 * @brief W_Prod Width to store a unormalised full product
		 */
		static constexpr unsigned int W_Prod = WE_Prod + WFF_Prod + 1;

		/**
		 * @brief PROD_MIN_LOWBIT_EXP weight of the low bit of the smallest non null product
		 */
		static constexpr int PROD_MIN_LOWBIT_EXP = 2*MIN_LOWBIT_WEIGHT;

		/**
		 * @brief PROD_MAX_EXP Maximum exponent of high significand bit of maximum normal value product
		 */
		static constexpr int PROD_MAX_EXP = MAX_NORMAL_UNBIASED_EXP * 2 + 1;

		/**
		 * @brief ACC_SIZE exact accumulator minimal width, including sign bit
		 */
		static constexpr unsigned int ACC_SIZE = static_cast<unsigned int>(PROD_MAX_EXP - PROD_MIN_LOWBIT_EXP + 1) + 1;

		/**
		 * @brief SMALLEST_SN_ACC_OFFSET index of the smallest representable subnormal weight in the accumulator
		 */
		static constexpr unsigned int MIN_SN_ACC_OFFSET = static_cast<unsigned int>(MIN_LOWBIT_WEIGHT - PROD_MIN_LOWBIT_EXP);

		/**
		 * @brief MAX_N_ACC_OFFSET index of the maximum representable weight in the accumulator
		 */
		static constexpr unsigned int MAX_N_ACC_OFFSET = static_cast<unsigned int>(MAX_NORMAL_UNBIASED_EXP - PROD_MIN_LOWBIT_EXP);

		/**
		 * @brief MIN_N_ACC_OFFSET index of the smallest normal bit in the accumulator
		 */
		static constexpr unsigned int MIN_N_ACC_OFFSET = MIN_SN_ACC_OFFSET + WF;
};

template<>
struct StandardIEEEDim<16> : public IEEEDim<5, 10> {};

template <>
struct StandardIEEEDim<32> : public IEEEDim<8, 23> {};

template<>
struct StandardIEEEDim<64> : public IEEEDim<11, 52> {};

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

			IEEENumber(Wrapper<WE + WF + 1, false> const & val = {0}):Wrapper<WE+WF+1, false>{val}{}

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
