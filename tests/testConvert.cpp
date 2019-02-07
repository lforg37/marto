#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE PositToValueTest

#include <cmath>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include "posit_decoder.hpp"
#include "posit_dim.hpp"

using namespace std;

int twoCompClean(int input, int nb_clear)
{
	int mask = -1;	
	int inverted = input ^ mask;
	int neg_val_dirty = inverted + 1; 
	int cleaned_val = neg_val_dirty & ((1 << nb_clear) - 1);
	return cleaned_val;
}

BOOST_AUTO_TEST_CASE(PositToValueTestPosit16) 
{
	for (size_t i = 0 ; i < (1<<16) - 1 ; ++i) {
		PositEncoding<16> encoding(i);
		auto decoded = posit_decoder(encoding);
		int value = decoded.getSignificand().to_int();
		if (decoded.getSignBit() == 1) {
			value = twoCompClean(value, PositDim<16>::WF + 1);
		}
		
		double floatingPointVal = value;

		int exp = decoded.getExp();
		if (exp >= (1 << (PositDim<16>::WE - 1))) { // exp should be negated
			exp = twoCompClean(exp, PositDim<16>::WE);
		}
		double factor = pow(2.0, exp);
		value *= factor;
		//TODO Compare with the corresponding posit N

	}
	BOOST_REQUIRE_MESSAGE(false, "test!!!");
}

