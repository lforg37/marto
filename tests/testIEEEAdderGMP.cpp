#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestIEEEAdder

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <gmpxx.h>

#include "hint.hpp"

#include "ieeefloats/ieee_adder.hpp"

#include "numeric_formats/ieee_small.hpp"

using namespace std;

BOOST_AUTO_TEST_CASE(TestIEEAddWE3WF4GMP)
{
	constexpr unsigned int WF = 4;
	constexpr unsigned int WE = 3;

	uint32_t nbRepr = uint64_t{1} << (WF+WE+1);

	libnumform::SmallIEEENumber<WE, WF> op1{0}, op2{0};


//	auto filepath = TESTDATA_ROOT"/ieeeadder/test_WE3_WF4.input";
//	ifstream testfile;
//	testfile.open(filepath);
//	string i0;
//	string i1;
//	string res;
//	BOOST_REQUIRE_MESSAGE(!testfile.fail(), "The test file " << filepath << " cannot be opened.");
	for (mpz_class op1_repr = 0 ; op1_repr < nbRepr ; ++op1_repr ) {
		op1 = static_cast<uint32_t>(op1_repr.get_ui());
		hint::IEEENumber<WE, WF, hint::GMPWrapper> op1_hint{op1_repr};
		for (mpz_class op2_repr = 0 ; op2_repr < nbRepr ; ++op2_repr ) {
			//cerr << "Iter " << op1_repr.get_str(2)  << ", " << op2_repr.get_str(2) << endl;
			op2 = static_cast<uint32_t>(op2_repr.get_ui());
			hint::IEEENumber<WE, WF, hint::GMPWrapper> op2_hint{op2_repr};

			auto adder_res = ieee_add_sub_impl(op1_hint, op2_hint);
			auto adder_res_gmp = adder_res.unravel();

			auto test_res = op1 + op2;
			if (test_res.isNaN()) {
				BOOST_REQUIRE_MESSAGE(adder_res.isNaN().template isSet<0>(), "Result should be NaN");
			} else {
				auto res_repr = mpz_class{test_res.getRepr()};
				BOOST_REQUIRE_MESSAGE(adder_res_gmp == res_repr, "Error for iteration " <<
								  op1_repr << ", " << op2_repr <<
								  "\nHint obtained result is\n" << adder_res_gmp.get_str(2) <<
								  "\nShould be :\n" << res_repr.get_str(2));
			}
		}
	}
}//*/
