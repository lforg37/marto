#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE PositToValueTest

#include <array>
#include <cmath>
#include <iostream>
#include <limits>

#include <boost/test/unit_test.hpp>
#ifdef SOFTPOSIT
#include "softposit.h"
#endif

#include "add_sub_quire.hpp"
#include "lzoc_shifter.hpp"
#include "posit_decoder.hpp"
#include "posit_dim.hpp"
#include "posit_encoder.hpp"
#include "posit_mul.hpp"
#include "quire_to_posit.hpp"
#include "value_prod_conversions.hpp"

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

BOOST_AUTO_TEST_CASE(PositValueToProd)
{
	uint16_t value = 0;
	PositValue<16> one(0, PositDim<16>::EXP_BIAS, 0, 1, 0);
	do {
		PositEncoding<16> posit_encoding{value};
		auto posit_val = posit_decoder<16>(posit_encoding);

		auto posit_prod_direct = PositValue_to_PositProd(posit_val);
		auto posit_prod_by_one = posit_mul<16>(posit_val, one);

		if (posit_val.getIsNaR() == 1) {
			BOOST_REQUIRE_MESSAGE(posit_prod_by_one.getIsNaR() == 1, "Prod by one should be NAR");
			BOOST_REQUIRE_MESSAGE(posit_prod_direct.getIsNaR() == 1, "Direct conversion should be NAR");
		} else {
			if (posit_prod_by_one != posit_prod_direct) {
				cerr << "=== posit_prod_by_one ===" << endl;
				posit_prod_by_one.printContent();
				cerr << "=== posit_prod_direct ===" << endl;
				posit_prod_direct.printContent();
			}
			BOOST_REQUIRE_MESSAGE(posit_prod_by_one == posit_prod_direct, "Error for conversion with value " << value);
		}
		value += 1;
	} while (value != 0);
}

BOOST_AUTO_TEST_CASE(PositValueToProdToValue)
{
	uint16_t value = 0;
	do {
		PositEncoding<16> current = value;
		auto decoded = posit_decoder(current);
		auto prod = PositValue_to_PositProd(decoded);
		auto casted_val = PositProd_to_PositValue(prod);
		auto reencoding = posit_encoder(casted_val);

		BOOST_REQUIRE_MESSAGE(reencoding == current, "Error for conversion with value " << value);		
		value += 1;
	} while (value != 0);
}

BOOST_AUTO_TEST_CASE(TestOppositeProd)
{
	PositValue<16> minus_one(
			0,
			28,
			1,
			0,
			0
		);
 	for (int16_t i = 0 ; i >= 0 ; ++i) {
		PositEncoding<16> enc{i};
		PositEncoding<16> opposite{i * -1};
		auto value = posit_decoder(enc);
		auto neg_prod = posit_mul(value, minus_one);
		auto op_value = posit_decoder(opposite);
		auto convert_op = PositValue_to_PositProd(op_value);
		if (convert_op != neg_prod) {
			cerr << "=== convert_op ===" << endl;
			convert_op.printContent();
			cerr << "=== neg_prod ===" << endl;
			neg_prod.printContent();
		}
		BOOST_REQUIRE_MESSAGE(
				convert_op == neg_prod,
				"Error with encoding " << i
			);
	}
}


BOOST_AUTO_TEST_CASE(TestZeroExpZero) 
{
	PositEncoding<16> val{0};
	auto decoded = posit_decoder(val);
	BOOST_REQUIRE_MESSAGE(decoded.getExp() == 0, 
			"Decoded biased exp of 0 should be zero after decoding"
		);
}

BOOST_AUTO_TEST_CASE(TestQuireConvertBack)
{
	Quire<16> quire{0};

	for(uint32_t value = 0; value < (1<<16); value++) {
		auto valueEncoding = PositEncoding<16> (value);
		auto decoded = posit_decoder(valueEncoding);
		auto prod = PositValue_to_PositProd(decoded);
		cerr << "Prod : " << prod << endl;
		Quire<16> quireConvert = add_sub_quire(quire, prod, 0);
		cerr << "Quire convert : " << quireConvert << endl;
		auto back_convert = quire_to_posit(quireConvert);
		back_convert.printContent();
		
		BOOST_REQUIRE_MESSAGE(back_convert == decoded,
				"Error for posit with encoding " << value
			);
	}
}

#ifdef SOFTPOSIT
BOOST_AUTO_TEST_CASE(TestMinMax)
{
	array<uint16_t, 4> encodings = {
		1, // MinPos
		(1<<15) - 1, //MaxPos
		numeric_limits<uint16_t>::max(), //MinNeg
		(1<<15) + 1 //MaxNeg
	};

	array<PositValue<16>, 4> values = {
		PositValue<16>::getMinPos(),
		PositValue<16>::getMaxPos(),
		PositValue<16>::getMinNeg(),
		PositValue<16>::getMaxNeg()
	};

	array<string, 4> names = {
		"MinPos", 
		"MaxPos",
		"MinNeg",
		"MaxNeg"
	};

	for (size_t i = 0 ; i < 4 ; ++i) {
		uint16_t enc = encodings[i];
		PositValue<16> & valStat = values[i];
		string& name = names[i];
		auto decoded = posit_decoder(PositEncoding<16>{enc});
		BOOST_REQUIRE_MESSAGE(decoded == valStat, 
				"Error in " << name << " generated value should be identical to decoded one");
	}

}

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
	// Should be tested on 32 bits as well when possible
	uint16_t i = 0;
	do {
		PositEncoding<16> encoding(i);
		auto decoded = posit_decoder(encoding);
		auto encoded = posit_encoder(decoded);

		BOOST_REQUIRE_MESSAGE(encoding == encoded,
				"Error : " << i << " is not encoded correctly");
		i += 1;
	} while (i != 0); 
}
