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

	using dim1 = FPDim<WE1, WF1>;
	using dim2 = FPDim<WE2, WF2>;

	using fpnum1 = FPNumber<dim1, VivadoWrapper>;

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
