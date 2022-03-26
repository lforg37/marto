#include <boost/test/unit_test.hpp>

namespace utf = boost::unit_test;

#ifdef MPFR
#include <mpfr.h>
#endif

#include <hint.hpp>

using hint::VivadoWrapper;

#include "floatingpoint/fp_number.hpp"


BOOST_AUTO_TEST_CASE(TestRounderLargerformat)
{
	constexpr vec_width WE1 = 4;
	constexpr vec_width WF1 = 7;

	constexpr vec_width WE2 = 5;
	constexpr vec_width WF2 = 9;

	using dim1 = FixedFormat<WE1, WF1>;
	using dim2 = FixedFormat<WE2, WF2>;

	using fpnum1 = FixedNumber<dim1, VivadoWrapper>;

	using rounder = Rounder<dim2, dim1, true>;

	for(uint64_t frac1 = 0 ; frac1 < (1<<WF1) ; ++frac1) {
		for (int64_t exp1 = dim1::MIN_EXP ; exp1 <= dim1::MAX_EXP ; exp1++) {
			for (uint32_t sign : {0, 1}){
				fpnum1 val{{frac1}, {exp1}, {sign}, {0}, {0}, {0}};
				auto rounded = rounder::compute(val);
				BOOST_REQUIRE_EQUAL(exp1, rounded.getExponent().unravel());
				BOOST_REQUIRE_EQUAL(frac1 << (WF2 - WF1), rounded.getFraction().unravel());
				BOOST_REQUIRE_EQUAL(sign, rounded.getSign().unravel());
				BOOST_REQUIRE_EQUAL(rounded.isInf().unravel(), 0);
				BOOST_REQUIRE_EQUAL(rounded.isNaN().unravel(), 0);
				BOOST_REQUIRE_EQUAL(rounded.isZero().unravel(), 0);
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(TestRounderSmallerformat)
{
	constexpr vec_width WE1 = 4;
	constexpr vec_width WF1 = 7;

	constexpr vec_width WE2 = 5;
	constexpr vec_width WF2 = 9;

	using dim1 = FixedFormat<WE1, WF1>;
	using dim2 = FixedFormat<WE2, WF2>;

	using fpnum = FixedNumber<dim2, VivadoWrapper>;

	using rounder = Rounder<dim1, dim2>;

	constexpr uint64_t overflow_mask{(1 << WF1)};

	bool all_ok = true;

	for(uint64_t frac = 0 ; frac < (1<<WF2) ; ++frac) {
		uint64_t truncated_frac = frac >> (WF2 - WF1 - 1);
		uint64_t r_frac = truncated_frac + 1;
		uint64_t expected_frac_wo = r_frac >> 1;
		bool overflow_frac = (expected_frac_wo & overflow_mask);
		uint64_t expected_frac_val = expected_frac_wo & (overflow_mask - 1);

		for (int64_t exp = dim2::MIN_EXP ; exp <= dim2::MAX_EXP ; exp++) {
			bool overflow_limit = (exp == dim1::MIN_EXP - 1);
			bool overflow = overflow_frac || overflow_limit;
			uint64_t expected_frac = (overflow_limit) ? 0 : expected_frac_val;
			int64_t expected_exp = exp  + ((overflow) ? 1 : 0);
			for (uint32_t sign : {0, 1}){
				fpnum val{{frac}, {exp}, {sign}, {0}, {0}, {0}};
				auto rounded = rounder::compute(val);
				auto rexp = rounded.getExponent().unravel();
				auto rfrac = rounded.getFraction().unravel();
				auto rs = rounded.getSign().unravel();
				auto isInf = rounded.isInf().unravel();
				auto isNan = rounded.isNaN().unravel();
				auto isZero = rounded.isZero().unravel();
				all_ok = all_ok and (rs == sign) and (isNan == 0);
				if  ((expected_exp <= dim1::MAX_STORABLE_EXP && expected_exp >= dim1::MIN_STORABLE_EXP) || overflow_limit) {
					all_ok = all_ok and (expected_frac == rfrac);
					all_ok = all_ok and (expected_exp == rexp);
					all_ok = all_ok and (isInf == 0);
					all_ok = all_ok and (isZero == 0);
				} else if (expected_exp > dim1::MAX_STORABLE_EXP) {
					all_ok = all_ok and (isInf == 1) and (isZero == 0);
				} else {
					all_ok = all_ok and (isInf == 0) and (isZero == 1);
				}

				if (not(all_ok)) {
					auto prettyprint = [](string const & title, uint64_t sign, int64_t exp, uint64_t frac) {
						cerr << title << ": <" << sign << ", " << exp << ", " << frac << ">\n";
					};
					prettyprint("InitValue", sign, exp, frac);
					prettyprint("Expecting", sign, expected_exp, expected_frac);
					prettyprint("Got", rs, rexp, rfrac);
					BOOST_REQUIRE(false);
				}
			}
		}
	}
}
