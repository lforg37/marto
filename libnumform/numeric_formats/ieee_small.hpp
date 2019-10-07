#ifndef IEEE_SMALL_HPP
#define IEEE_SMALL_HPP

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
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
			SmallIEEENumber(double upper_half, double lower_half):SmallIEEENumber{}
			{
				if (std::isnan(upper_half)) {
					_repr = _NaN;
					return;
				}

				if (std::isinf(upper_half)) {
					_repr = _exp_mask;
					if (upper_half < 0)
						_repr |= _sign_mask;
					return;
				}

				if (upper_half == .0) {
					_repr = 0;
					if (std::signbit(upper_half))
						_repr |= _sign_mask;
					return;
				}
				int exp_up, exp_low;
				bool neg_up, neg_low;
				uint32_t signif_up, signif_low;
				neg_up = upper_half < 0;
				if  (neg_up) {
					upper_half *= -1.0;
				}

				neg_low = lower_half < 0;
				if (neg_low) {
					lower_half *= -1;
				}

				bool same_sign = (neg_low == neg_up);
				constexpr double EXP_SCALE_DOUBLE = static_cast<double>(uint64_t{1} << 53);

				uint64_t up_significand = static_cast<uint64_t>(frexp(upper_half, &exp_up) * EXP_SCALE_DOUBLE);
				exp_up -= 1;

				if (exp_up > _emax) {
					_repr = (static_cast<uint32_t>(neg_up) << (WF + WE)) | _exp_mask;
					return;
				}

				constexpr uint8_t left_out_bit_count = 52 - WF;
				constexpr uint64_t round_mask = uint64_t{1} << (left_out_bit_count - 1);
				constexpr uint64_t sticky_mask = round_mask - 1;
				bool round_bit = round_mask & up_significand;
				bool sticky_up = sticky_mask & up_significand;
				bool sticky_down = lower_half > .0;
				uint8_t extra_info_vec = (round_bit << 2) | (sticky_up << 1);
				if (same_sign) {
					extra_info_vec |= sticky_down;
				} else {
					extra_info_vec -= sticky_down;
				}

				up_significand >>= left_out_bit_count;

				bool sticky = extra_info_vec & 3;
				round_bit = extra_info_vec & 4;
				bool borrow = extra_info_vec & 8;

				if (borrow) {
					up_significand -= 1;
				}

				constexpr uint32_t underflow_mask = uint32_t{1} << WF;

				if (not (underflow_mask & up_significand)) {
					exp_up -= 1;
				}

				if (exp_up < _emin) {
					int shifts = _emin - exp_up;
					shifts = (shifts > WF + 2) ? WF + 2 : shifts;
					uint64_t round_mask = (uint64_t{1} << (shifts-1));
					sticky = sticky or round_bit or ((round_mask - 1) & up_significand);
					round_bit = up_significand & round_mask;
					exp_up = _emin - 1;
					up_significand >>= shifts;
				}

				exp_up += _bias;
				uint32_t value = (static_cast<uint32_t>(up_significand) & _fraction_mask) |
						(exp_up << WF) |
						(static_cast<uint32_t>(neg_up) << (WF + WE));

				if (round_bit and (sticky or (value & 1))) {
					value += 1;
				}
				_repr = value;
			}

			double getBinary64() {
				if ((_repr & _exp_mask) == _exp_mask) {
					if (_repr & _fraction_mask) {
						return std::numeric_limits<double>::quiet_NaN();
					} else {
						constexpr double inf = std::numeric_limits<double>::infinity();
						constexpr double m_inf = -inf;
						return (_repr & _sign_mask) ? m_inf : inf;
					}
				}
				int32_t exp = getExp();
				exp += (exp == 0);
				exp -=  _bias + WF;
				int frac = getExplicitFrac();
				double res = ldexp(frac, exp);
				if (_repr & _sign_mask)
					res *= -1.0;
				return res;
			}

			SmallIEEENumber(uint32_t repr):SmallIEEENumber{}
			{
				_repr = repr;
				if (DynamicValueCheck) {
					if ((repr & (~_format_mask)) != 0) {
						std::cerr << "Uncorrect initialisation value " << repr << std::endl;
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
						std::cerr << "Uncorrect initialisation value " << repr << std::endl;
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
			static constexpr int32_t _emin = -1 * (int32_t{1} << (WE - 1)) + 2;
			static constexpr int32_t _emax = (uint32_t{1} << (WE - 1)) - 1;

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
		auto op1_val = op1.getBinary64();
		auto op2_val = op2.getBinary64();

		double prod = op1_val * op2_val;
		double err = std::fma(op1_val, op2_val, -prod);

		return {prod, err};
	}

	template<unsigned int WE, unsigned int WF, bool DynamicValueCheck>
	SmallIEEENumber<WE, WF, DynamicValueCheck> operator+(
			SmallIEEENumber<WE, WF, DynamicValueCheck> op1,
			SmallIEEENumber<WE, WF, DynamicValueCheck> op2
		)
	{
		using SIEEE = SmallIEEENumber<WE, WF, DynamicValueCheck>;
		double op1_val = op1.getBinary64();
		double op2_val = op2.getBinary64();

		double s = op1_val + op2_val;
		double op1_approx = s - op2_val;
		double op2_approx = s - op1_approx;
		double d_op1 = op1_val - op1_approx;
		double d_op2 = op2_val - op2_approx;
		double err = d_op1 + d_op2;
		return  {s, err};
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
