#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE PositToValueTest

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/test/unit_test.hpp>
#ifdef SOFTPOSIT
#include "softposit.h"
#endif

#include "posit_decoder.hpp"
#include "posit_encoder.hpp"
#include "posit_dim.hpp"
#include "lzoc_shifter.hpp"

using namespace std;

int twoCompClean(int input, int nb_clear)
{
	int mask = -1;	
	int inverted = input ^ mask;
	int neg_val_dirty = inverted + 1; 
	int cleaned_val = neg_val_dirty & ((1 << nb_clear) - 1);
	return -1 * cleaned_val;
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

#ifdef SOFTPOSIT
BOOST_AUTO_TEST_CASE(PositToValueTestPosit16) 
{
	uint16_t i = 0;
	do {
		//cout << "Testing " << i << endl;
		PositEncoding<16> encoding(i);
		auto decoded = posit_decoder(encoding);
		ap_uint<14> complete_sig = decoded.getSignBit().concat(decoded.getSignificand());
		int value = complete_sig.to_int();
		if (decoded.getSignBit() == 1) {
			value = twoCompClean(value, PositDim<16>::WF + 2);
		}
		
		double floatingPointVal = value;

		int exp = decoded.getExp() - PositDim<16>::EXP_BIAS;
		double factor = pow(2.0, (exp-12));
		floatingPointVal *= factor;
		//cout << "Decoded value : " << value << endl;
		
		if (decoded.getIsNaR() == 1) {
			floatingPointVal = numeric_limits<double>::infinity();
		}

		posit16_t posit_val = castP16(i);
		double soft_posit_val = convertP16ToDouble(posit_val);
		//cout << "Soft Posit decoded value : " << soft_posit_val << endl << endl;
		BOOST_REQUIRE_MESSAGE(
				soft_posit_val == floatingPointVal,
			   "Error in conversion : decoding of value " << i << 
			   " gives " << floatingPointVal << " instead of " << soft_posit_val );
		i += 1;
	}while(i != 0);
}
#endif


BOOST_AUTO_TEST_CASE(PositEncoder)
{
	uint16_t i = 0;
	do {
		PositEncoding<16> encoding(i);
		auto decoded = posit_decoder(encoding);
		auto encoded = posit_encoder(decoded);

		BOOST_REQUIRE_MESSAGE(encoding == encoded,
				"Error : " << i << " gave is not encoded correctly");
		i += 1;
	} while (i != 0); 
}
