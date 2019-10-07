#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestLibNumForm

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "numeric_formats/ieee_small.hpp"

using namespace std;

BOOST_AUTO_TEST_CASE(testConvertBack)
{
	constexpr unsigned int WE = 3;
	constexpr unsigned int WF = 4;

	using SIEEE = libnumform::SmallIEEENumber<WE, WF>;

	for (uint32_t t_idx = 0 ; t_idx < (1 << (WE + WF + 1)) ; ++t_idx) {
		SIEEE val{t_idx};
		double test = val.getBinary64();
		SIEEE val_comp(test, 0.);
		if (val.isNaN()) {
			BOOST_REQUIRE_MESSAGE(val_comp.isNaN(), "Error at iteration " << t_idx <<
								  ". Expecting NaN, got something different");
		} else {
			BOOST_REQUIRE_MESSAGE(val_comp.getRepr() == val.getRepr(),
				"Error for iteration " << t_idx << ". Expecting " << val.getRepr() <<
				" got " << val_comp.getRepr());

		}
	}
}
