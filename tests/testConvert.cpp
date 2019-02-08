#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE PositToValueTest

#include <cmath>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include "softposit.h"

#include "posit_decoder.hpp"
#include "posit_dim.hpp"
#include "lzoc_shifter.hpp"

using namespace std;

int twoCompClean(int input, int nb_clear)
{
	int mask = -1;	
	int inverted = input ^ mask;
	int neg_val_dirty = inverted + 1; 
	int cleaned_val = neg_val_dirty & ((1 << nb_clear) - 1);
	return cleaned_val;
}

BOOST_AUTO_TEST_CASE(LZOCShiftTest)
{
	uint16_t i = 0;
	do {
		ap_uint<1<<4> entry = i;
		auto ret = lzoc_shifter<4>(entry, ap_uint<1>(0));	
		auto shift = ret.range(19, 16);
		auto val = ret.range(15, 0);
		unsigned int intval = val.to_uint();
		unsigned int shiftint = shift.to_uint();
		BOOST_REQUIRE_MESSAGE((intval >> shiftint) == (unsigned int) i,
				"Error : " << i << " gave " << ret);
		i += 1;
	} while (i != 0); 
}

BOOST_AUTO_TEST_CASE(PositToValueTestPosit16) 
{
	uint16_t i = 0;
	do {
		//cout << "Testing " << i << endl;
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
		//cout << "Decoded value : " << value << endl;

		posit16_t posit_val = castP16(i);
		double soft_posit_val = convertP16ToDouble(posit_val);
		//cout << "Soft Posit decoded value : " << soft_posit_val << endl << endl;
		BOOST_REQUIRE_MESSAGE(
				soft_posit_val == value,
			   "Error in conversion : decoding of value " << i << 
			   " gives " << value << " instead of " << soft_posit_val );
		i += 1;
	}while(i != 0);
}

