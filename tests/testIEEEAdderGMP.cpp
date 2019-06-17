#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestIEEEAdder

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <gmpxx.h>

#include "hint.hpp"

#include "ieeefloats/ieee_adder.hpp"

using namespace std;

BOOST_AUTO_TEST_CASE(TestIEEAddWE3WF4GMP)
{
	constexpr unsigned int WF = 4;
	constexpr unsigned int WE = 3;

	uint64_t nbTests = 1 << ((WF+WE+1) * 2);

	auto filepath = TESTDATA_ROOT"/ieeeadder/test_WE3_WF4.input";
	ifstream testfile;
	testfile.open(filepath);
	string i0;
	string i1;
	string res;
	BOOST_REQUIRE_MESSAGE(!testfile.fail(), "The test file " << filepath << " cannot be opened.");
	for (uint64_t curtest = 0; curtest < nbTests ; ++curtest) {
		testfile >> i0 >> i1 >> res >> res;
		mpz_class i0gmp{i0, 2};
		mpz_class i1gmp{i1, 2};
		mpz_class mpzres{res, 2};
		IEEENumber<WE, WF, hint::GMPWrapper> i0ieee{i0gmp};
		IEEENumber<WE, WF, hint::GMPWrapper> i1ieee{i1gmp};

		auto result = ieee_add_sub_impl(i0ieee, i1ieee);

		//cerr << "IN0 : " << endl << i0gmp.get_str(2) << endl;
		//cerr << "IN1 : " << endl << i1gmp.get_str(2) << endl;
		//cerr << "EXPECTED : " << endl << mpzres.get_str(2) << endl;
		//cerr << "GOT : " << endl << result.unravel().get_str(2) << endl;

		BOOST_REQUIRE_MESSAGE(result.unravel() == mpzres, "Error for iteration " << curtest);
	}
}
