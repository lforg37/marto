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

			us_wrapper<1> getSign() const{
				return basetype::template get<WE+WF>();
			}

			us_wrapper<WE> getExponent() const {
				return basetype::template slice<WF + WE - 1, WF>();
			}

			us_wrapper<WF> getFractionnalPart() const {
				return basetype::template slice<WF - 1, 0>();
			}

			us_wrapper<1> getLeadBitVal() const {
				return getExponent().or_reduction();
			}

			us_wrapper<1> isInfinity() const {
				return getExponent().and_reduction().bitwise_and(getFractionnalPart().or_reduction().invert());
			}

			us_wrapper<1> isNaN() const {
				return getExponent().and_reduction().bitwise_and(getFractionnalPart().or_reduction());
			}
	};

	/*
	template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
	class IEEEProduct: public Wrapper<WE + 1 + 2*WF + 2 + 1, false>
	{
		private:
			static constexpr unsigned int encoding_size = 1 + (WE + 1) + 2*(WF+1); // Sign, exp, explicited fraction
			typedef Wrapper<WE+WF+1, false> basetype;
			template<unsigned int W>
			using us_wrapper = Wrapper<W, false>;
	};*/

#endif // IEEETYPE_HPP
