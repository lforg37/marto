#ifndef IEEE_SMALL_HPP
#define IEEE_SMALL_HPP

#include <cstdint>
#include <type_traits>

namespace libnumform {

	template<typename Type, unsigned int NbShift, unsigned int MaxSize, typename Activate=void>
	struct _CLZOnePropagator;

	template<typename Type, unsigned int NbShift, unsigned int MaxSize>
	struct _CLZOnePropagator<Type, NbShift, MaxSize, typename std::enable_if<(NbShift<=MaxSize)>::type>
	{
		static inline void propagate(Type& value)
		{
			value |= value >> NbShift;
			constexpr unsigned int next_level_shift = NbShift * 2;
			_CLZOnePropagator<Type, next_level_shift, MaxSize>::propagate(value);
		}
	};

	template<typename Type, unsigned int NbShift, unsigned int MaxSize>
	struct _CLZOnePropagator<Type, NbShift, MaxSize, typename std::enable_if<(NbShift>MaxSize)>::type>
	{
		static inline void propagate(Type&)
		{
			return;
		}
	};

	template<typename Type, Type Mask, uint8_t Cost, typename Activate = void>
	struct _CLZTestDec;

	template<typename Type, Type Mask, uint8_t Cost>
	struct _CLZTestDec<Type, Mask, Cost, typename std::enable_if<(Cost >= 1)>::type>
	{
		static inline void test_dec_count(Type const & value, uint8_t & count)
		{
			if (value & Mask) count -= Cost;
			constexpr uint8_t next_level_cost = Cost >> 1;
			constexpr Type next_level_mask = Mask ^ (Mask >> next_level_cost);
			_CLZTestDec<Type, next_level_mask, next_level_cost>::test_dec_count(value, count);
		}
	};

	template<typename Type, Type Mask, uint8_t Cost>
	struct _CLZTestDec<Type, Mask, Cost, typename std::enable_if<(Cost == 0)>::type>
	{
		static inline void test_dec_count(Type const &, uint8_t &)
		{
			return;
		}
	};

	template<typename Type>
	inline uint8_t clz(Type val){
		static_assert(std::is_unsigned<Type>::value and std::is_integral<Type>::value, "clz should be called on an unsigned integer");
		uint8_t count = sizeof(Type) << 3;
		if (val == 0)
			return count;
		count -= 1;
		_CLZOnePropagator<Type, 1, sizeof(Type) * 4>::propagate(val);
		val >>= 1;
		val += 1;
		constexpr Type mask = ~Type{0} ^ ((~Type{0}) >> (sizeof(Type) << 2));
		_CLZTestDec<Type, mask, sizeof(Type) * 4>::test_dec_count(val, count);

		return count;
	}

	template<unsigned int WE, unsigned int WF, bool DynamicValueCheck>
	class SmallIEEENumber;

	template<unsigned int WE, unsigned int WF, bool DynamicValueCheck>
	SmallIEEENumber<WE, WF, DynamicValueCheck> operator*(
			SmallIEEENumber<WE, WF, DynamicValueCheck> op1,
			SmallIEEENumber<WE, WF, DynamicValueCheck> op2
		);

	template<unsigned int WE, unsigned int WF, bool DynamicValueCheck>
	SmallIEEENumber<WE, WF, DynamicValueCheck> operator+(
			SmallIEEENumber<WE, WF, DynamicValueCheck> op1,
			SmallIEEENumber<WE, WF, DynamicValueCheck> op2
		);

	template<unsigned int WE, unsigned int WF, bool DynamicValueCheck>
	SmallIEEENumber<WE, WF, DynamicValueCheck> operator-(
			SmallIEEENumber<WE, WF, DynamicValueCheck> op1,
			SmallIEEENumber<WE, WF, DynamicValueCheck> op2
		);

	template<unsigned int WE, unsigned int WF, bool DynamicValueCheck = true>
	class SmallIEEENumber
	{
		private:
			SmallIEEENumber()
			{
				static_assert(WE >= 3, "SmallIEEENumber requires a exponent field of at least 3 bits");
				static_assert(WE + WF <= 31, "The type cannot handle floating point formats on more than 32 bits");
				static_assert(WF > 0, "Fraction field size should be at least 1");
			}

			inline uint32_t getExplicitFrac() const
			{
				return (_repr & _fraction_mask) | (isNullExp() ? 0 : _implicit_bit);
			}

			inline uint32_t getSign() const
			{
				return _repr & _sign_mask;
			}

			inline uint32_t getExp() const
			{
				return (_repr & _exp_mask) >> WF;
			}

			inline SmallIEEENumber getOpposite() const
			{
				return _repr ^ _sign_mask;
			}



		public:
			SmallIEEENumber(uint32_t repr):SmallIEEENumber{}
			{
				_repr = repr;
				if (DynamicValueCheck) {
					if ((repr & (~_format_mask)) != 0) {
						throw "Error, trying to initialise IEEENumber with a value greater than the format width";
					}
				}
			}

			inline bool isSubnormal() const
			{
				return ((_repr & _exp_mask) == 0 ) and ((_repr & _fraction_mask) != 0);
			}

			inline bool isZero() const
			{
				constexpr uint32_t inv_sign_mask = ~ _sign_mask;
				return (_repr & inv_sign_mask) == 0;
			}

			inline bool isNullExp() const
			{
				return (_repr & _exp_mask) == 0;
			}

			inline bool isInfinite() const
			{
				constexpr uint32_t inv_sign_mask = ~ _sign_mask;
				return (_repr & inv_sign_mask) == _exp_mask;
			}

			inline bool isNaN() const
			{
				return ((_repr & _exp_mask) == _exp_mask) and ((_repr & _fraction_mask) != 0);
			}

			SmallIEEENumber & operator=(uint32_t repr)
			{
				if (DynamicValueCheck) {
					if ((repr & (~_format_mask)) != 0) {
						throw "Error, trying to initialise IEEENumber with a value greater than the format width";
					}
				}
				_repr = repr;

				return *this;
			}

			inline uint32_t getRepr() const
			{
				return _repr;
			}

		private:
			uint32_t _repr;
			static constexpr unsigned int _bias = (1 << (WE - 1)) - 1;
			static constexpr uint32_t _exception_exp = (uint32_t{1} << WE) - 1;
			static constexpr unsigned int _format_size = WE + WF + 1;
			static constexpr uint32_t _exp_mask = ((1 << WE) - 1) << WF;
			static constexpr uint32_t _fraction_mask = (1 << WF) - 1;
			static constexpr uint32_t _sign_mask = 1 << (WE + WF);
			static constexpr uint32_t _implicit_bit = 1 << WF;
			static constexpr uint32_t _NaN = _exp_mask | 1;
			static constexpr uint32_t _format_mask = _sign_mask | _fraction_mask | _exp_mask;

			friend SmallIEEENumber operator*<WE, WF, DynamicValueCheck>(SmallIEEENumber op1, SmallIEEENumber op2);
			friend SmallIEEENumber operator+<WE, WF, DynamicValueCheck>(SmallIEEENumber op1, SmallIEEENumber op2);
			friend SmallIEEENumber operator-<WE, WF, DynamicValueCheck>(SmallIEEENumber op1, SmallIEEENumber op2);
	};

	template<unsigned int WE, unsigned int WF, bool DynamicValueCheck>
	SmallIEEENumber<WE, WF, DynamicValueCheck> operator*(
			SmallIEEENumber<WE, WF, DynamicValueCheck> op1,
			SmallIEEENumber<WE, WF, DynamicValueCheck> op2
		)
	{
		using SIEEE = SmallIEEENumber<WE, WF, DynamicValueCheck>;
		uint32_t sign = op1.getSign() ^ op2.getSign();
		bool op1_nan = op1.isNaN();
		bool op2_nan = op2.isNaN();
		bool op1_infinite = op1.isInfinite();
		bool op2_infinite = op2.isInfinite();
		bool op1_zero = op1.zero();
		bool op2_zero = op2.isZero();
		bool isResultNan = op1_nan or op2_nan or ((op1_infinite xor op2_infinite) and (op1_zero xor op2_zero));
		if (isResultNan)
			return {SIEEE::_NaN | sign};
		if (op1_zero or op2_zero)
			return {sign};
		if (op1_infinite or op2_infinite)
			return {SIEEE::_exp_mask | sign};
		uint64_t op1_significand = static_cast<uint64_t>(op1.getExplicitFrac());
		uint64_t op2_significand = static_cast<uint64_t>(op2.getExplicitFrac());
		uint64_t product = op1_significand * op2_significand;

		constexpr uint64_t prod_overflow_mask = uint64_t{1} << (2*WF + 1);
		bool prod_overflow = product & prod_overflow_mask;

		bool sticky = false;
		if (prod_overflow) {
			sticky = product & 1;
			product >>= 1;
		}

		int32_t exp = static_cast<int32_t>(op1.getExp()) +
				static_cast<int32_t>(op2.getExp()) -
				static_cast<int32_t>(SIEEE::_bias) + prod_overflow;

		constexpr int32_t subnormal_emin = -WF;
		constexpr uint64_t normal_mask = prod_overflow_mask >> 1;
		constexpr uint64_t normal_frac_mask = normal_mask - 1;
		constexpr uint64_t fraction_mask = ((uint64_t{1} << WF) - 1) << WF;

		bool round_bit;
		uint32_t value;

		if (exp >= 0) {
			constexpr uint8_t known_zeros = 64 - (2 * WF + 1);
			uint8_t prod_clz = clz(product) - known_zeros;

			uint8_t norm_budget = static_cast<uint8_t>((exp > prod_clz) ? prod_clz : exp);
			exp -= norm_budget;
			product <<= norm_budget;
			uint64_t fraction = (product & fraction_mask) >> WF;
			constexpr uint64_t round_mask = 1 << (WF-1);
			round_bit = round_mask & product;
			constexpr uint64_t sticky_mask = round_mask - 1;
			bool sticky = sticky or (product & sticky_mask);
			value  = static_cast<uint32_t>(fraction) | sign | (static_cast<uint32_t>(exp) << WF);
		} else {
			if (exp < subnormal_emin) {
				bool round_up = (exp == subnormal_emin) and (product & normal_mask) and (product & normal_frac_mask);
				return {sign | uint32_t{round_up}};
			}
			int32_t shift_budget = -exp;
			uint64_t round_mask = uint64_t{1} << (shift_budget - 1);
			 round_bit = product & round_mask ;
			round_mask -= 1;
			sticky = sticky or (product & round_mask);
			product >>= shift_budget;
			uint64_t fraction = (product & fraction_mask) >> WF;
			value = sign | static_cast<uint32_t>(fraction);

		}
		bool odd = value & 1;
		if (round_bit and (sticky or odd))
			value += 1;
		return {value};
	}

	template<unsigned int WE, unsigned int WF, bool DynamicValueCheck>
	SmallIEEENumber<WE, WF, DynamicValueCheck> operator+(
			SmallIEEENumber<WE, WF, DynamicValueCheck> op1,
			SmallIEEENumber<WE, WF, DynamicValueCheck> op2
		)
	{
		using SIEEE = SmallIEEENumber<WE, WF, DynamicValueCheck>;
		auto op1_sign = op1.getSign();
		auto op2_sign = op2.getSign();
		bool opposite_sign = op1_sign ^ op2_sign;
		bool op1_nan = op1.isNaN();
		bool op2_nan = op2.isNaN();
		bool op1_infinite = op1.isInfinite();
		bool op2_infinite = op2.isInfinite();

		if (op1.isZero()) {
			if(op2.isZero()) {
				return {op1_sign & op2_sign};
			}
			return op2;
		}

		if (op2.isZero()) {
			return op1;
		}

		if (op1_nan or op2_nan or (op1_infinite and op2_infinite and opposite_sign)) {
			return {SIEEE::_NaN};
		}

		if ((op1.getRepr() ^ op2.getRepr()) == SIEEE::_sign_mask) {
			return {0};
		}

		if (op1_infinite) {
			return {op1_sign | SIEEE::_exp_mask};
		}

		if (op2_infinite) {
			return {op2_sign | SIEEE::_exp_mask};
		}

		auto exp1 = op1.getExp() ;
		auto exp2 = op2.getExp();

		auto frac1 = op1.getExplicitFrac();
		auto frac2 = op2.getExplicitFrac();

		if (exp2 > exp1 or ((exp1 == exp2) and frac2 > frac1)) {
			std::swap(op1, op2);
			std::swap(exp1, exp2);
			std::swap(op1_sign, op2_sign);
			std::swap(frac1, frac2);
		}

		auto diffexp = exp1  +  (exp1 == 0) - (exp2 + (exp2 == 0));
		if (diffexp > (WF+1)) {
			return op1;
		}

		bool round_bit = false;
		bool sticky = false;

		uint32_t value;

		if (not opposite_sign) {
			uint32_t sum_op2 = frac2 >> diffexp;
			uint32_t sum = frac1 + sum_op2;

			if (exp1 == 0) {
				return {sum | op1_sign};
			}
			constexpr uint32_t overflow_mask = uint32_t{1} << (WF + 1);
			bool overflow = overflow_mask & sum;
			bool guard = false;
			if (diffexp >= 1) {
				uint32_t guard_mask = uint32_t{1} << (diffexp - 1);
				guard = frac2 & guard_mask;
				sticky = frac2 & (guard_mask - 1);
			}
			if (overflow) {
				exp1 += 1;
				if (exp1 == SIEEE::_exception_exp) {
					return {op1_sign | SIEEE::_exp_mask};
				}
				round_bit = sum & 1;
				sticky = sticky or guard;
				sum >>= 1;
			} else {
				round_bit = guard;
			}
			value = op1_sign | (exp1 << WF) | (sum & SIEEE::_fraction_mask);
			if (round_bit and (sticky or (value & 1))) {
				value += 1;
			}
			return {value};
		} else {
			uint32_t sum = frac1;
			uint32_t sum_op2 = frac2 >> diffexp;
			sum -= sum_op2;

			if (diffexp > 1) { // No cancellation possible
				uint32_t guard_mask = uint32_t{1} << (diffexp - 1);
				bool guard = guard_mask & frac2;
				guard_mask >>= 1;
				bool guard2 = guard_mask & frac2;
				sticky = frac2 & (guard_mask - 1);

				int8_t val = (guard << 2) | (guard2 << 1) | sticky;

				bool remove_one = guard or sticky or guard2;
				if (remove_one) {
					sum -= 1;
					val = -val;
				}

				guard = val & 4;
				guard2 = val & 2;
				sticky = sticky & 1;

				constexpr uint32_t underflow_mask = uint32_t{1} << (WF);
				bool underflow = not(underflow_mask & sum);

				if (underflow) {
					exp1 -= 1; //Always possible as exp1 >= 2;
					sum <<= 1;
					sum |= guard;
					round_bit = guard2;
				} else {
					round_bit = guard;
					sticky = sticky or guard2;
				}

				auto frac_res = sum & SIEEE::_fraction_mask;
				value = op1_sign | (exp1 << WF) | frac_res;
				bool frac_is_zero = (frac_res == 0);

				bool odd = value & 1;

				if (round_bit & (sticky or odd)) {
					value += 1;
				}

				return {value};
			} else {
				constexpr uint32_t known_zeros = 32 - (1 + WF);
				round_bit = frac2 & (diffexp != 0);
				if (round_bit)
					sum -= 1;
				auto clz_sum = static_cast<uint32_t>(clz(sum));
				clz_sum -= known_zeros;
				auto shift_budget = (exp1 < clz_sum) ? exp1 : clz_sum;
				exp1 -= shift_budget;
				if ((exp1 == 0) and (shift_budget != 0)) shift_budget -= 1;
				if (shift_budget) {
					sum <<= shift_budget;
					sum |= round_bit << (shift_budget - 1);
					round_bit = false;
				}
				auto frac_res = sum & SIEEE::_fraction_mask;
				value = op1_sign | (exp1 << WF) | (sum & frac_res);

				if (round_bit and (sum & 1)) {
					value += 1;
				}
				return {value};
			}
		}

		/*
		uint64_t summand = static_cast<uint64_t>(frac1) << diffexp;
		uint64_t second_operand = static_cast<uint64_t>(frac2);

		bool round_bit;
		bool sticky;

		uint32_t value;

		if (opposite_sign) {
			uint64_t res = summand - second_operand;
			if (diffexp <= 1) {
				if ((diffexp == 0) and(frac1 == frac2))
					return 0;
				uint32_t known_zeros = 64 - (WF + diffexp + 1);
				uint32_t clz_frac = static_cast<uint32_t>(clz(res));
				clz_frac -= known_zeros;
				uint32_t shift_budget = ((exp1) < clz_frac) ? (exp1) : clz_frac;
				res <<= shift_budget;
				exp1 -= shift_budget;
			}


			if (diffexp == 0) { // No rounding bits
				return {op1_sign | (exp1 << WF) | (static_cast<uint32_t>(res) & SIEEE::_fraction_mask)};
			} else {
				uint64_t round_bit_mask = uint64_t{1} << (diffexp - 1);
				uint64_t sticky_mask = round_bit_mask - 1;
				round_bit = res & round_bit_mask;
				sticky = res & sticky_mask;
				value = op1_sign | (exp1 << WF) | (static_cast<uint32_t>(res >> diffexp) & SIEEE::_fraction_mask);
			}

		} else {
			summand += second_operand;
			if (exp1 == 0) {
				constexpr uint64_t res_mask = static_cast<uint64_t>(SIEEE::_format_mask ^ SIEEE::_sign_mask);
				return {op1_sign | static_cast<uint32_t>(res_mask & summand)};
			}
			uint64_t overflow_mask = (uint64_t{1} << (diffexp + WF + 1));
			bool overflow = overflow_mask & summand;
			if (not overflow) {
				summand <<= 1;
			} else {
				exp1 += 1;
			}
			uint64_t shifted_wf = static_cast<uint64_t>(SIEEE::_fraction_mask) << (diffexp + 1);
			uint64_t unrounded_frac = (shifted_wf & summand) >> (diffexp + 1);

			uint64_t round_bit_mask = uint64_t{1} << diffexp;
			uint64_t sticky_mask = round_bit_mask - 1;
			round_bit = round_bit_mask & summand;
			sticky = sticky_mask & summand;
			value = op1_sign | (exp1 << WF) | unrounded_frac;
		}

		bool odd = value & 1;

		if (round_bit and (sticky or odd)) {
			value += 1;
		}

		return {value};*/
	}

	template<unsigned int WE, unsigned int WF, bool DynamicValueCheck>
	SmallIEEENumber<WE, WF, DynamicValueCheck> operator-(
			SmallIEEENumber<WE, WF, DynamicValueCheck> op1,
			SmallIEEENumber<WE, WF, DynamicValueCheck> op2
		)
	{
		using SIEEE = SmallIEEENumber<WE, WF, DynamicValueCheck>;
		return op1 + SIEEE{op2.getRepr() ^ SIEEE::_sign_mask};
	}
}

#endif // IEEE_SMALL_HPP
