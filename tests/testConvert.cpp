#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE PositToValueTest

#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>

#include <boost/test/unit_test.hpp>
#ifdef SOFTPOSIT
#include "softposit.h"
#endif

#include "posit.hpp"

using namespace std;

int twoCompClean(int input, int nb_clear)
{
	int mask = -1;	
	int inverted = input ^ mask;
	int neg_val_dirty = inverted + 1; 
	int cleaned_val = neg_val_dirty & ((1 << nb_clear) - 1);
	return -1 * cleaned_val;
}

BOOST_AUTO_TEST_CASE(QuireBackCornerCases)
{
	//Positive underflow
    StandardQuire<16> quire{0};
    auto minpos = StandardPIF<16>::getMinPos();
	auto prod = posit_mul(minpos, minpos);
	auto quire_conv = add_sub_quire(quire, prod, 0);
	auto decoded = quire_to_posit(quire_conv);
	BOOST_REQUIRE_MESSAGE(decoded == minpos, 
			"Positive underflow does not returns minpos"
		);

	//Positive sticky
	prod = PositIF_to_PositProd(minpos);
	quire_conv = add_sub_quire(quire_conv, prod, 0);
	decoded = quire_to_posit(quire_conv);
    BOOST_REQUIRE_MESSAGE(decoded xor minpos == (1 << 19),
			"Error, sticky bit is not set for minpos * (1+minpos)"
		);

	//Positive overflow
    auto maxpos = StandardPIF<16>::getMaxPos();
	prod = posit_mul(maxpos, maxpos);
	quire_conv = add_sub_quire(quire, prod, 0);
	decoded = quire_to_posit(quire_conv);
	BOOST_REQUIRE_MESSAGE(decoded == maxpos,
			"Positive overflow doesn't return maxpos"
			);

	//Negative underflow
    auto minneg = StandardPIF<16>::getMinPos();
	prod = posit_mul(minpos, minneg);
	quire_conv = add_sub_quire(quire, prod, 0);
	decoded = quire_to_posit(quire_conv);
	BOOST_REQUIRE_MESSAGE(decoded == minneg,
			"Negative underflow is not mapped to minneg"
			);

	//Negative overflow
    auto maxneg = StandardPIF<16>::getMaxPos();
	prod = posit_mul(maxpos, maxneg);
	quire_conv = add_sub_quire(quire, prod, 0);
	decoded = quire_to_posit(quire_conv);
	BOOST_REQUIRE_MESSAGE(decoded == maxneg,
			"Negative overflow is not mapped to maxneg"
			);
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
    StandardPIF<16> one(0, StandardPositDim<16>::EXP_BIAS, 0, 1, 0);
	do {
        PositEncoding<16, 1> posit_encoding{value};
        auto posit_val = static_cast<PositIntermediateFormat<16, 1> >(posit_encoding);

        auto posit_prod_direct = static_cast<PositProd<16, 1> >(posit_val);
        auto posit_prod_by_one = posit_mul(posit_val, one);

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
        StandardPositEncoding<16> current{value};
        auto decoded = static_cast<StandardPIF<16> >(current);
        auto prod = static_cast<StandardPositProd<16> >(current);
		auto casted_val = PositProd_to_PositIF(prod);
        auto reencoding = static_cast<StandardPositEncoding<16> >(casted_val);

		if(reencoding != current){
			fprintf(stderr, "=== decoded ===\n");
			decoded.printContent();
			fprintf(stderr, "=== prod ===\n");
			prod.printContent();
			fprintf(stderr, "=== casted_val ===\n");
			casted_val.printContent();
            printApUint(reencoding);
            fprintf(stderr, "\n");
            printApUint(current);
           fprintf(stderr, "%d\n", value);
		}

		BOOST_REQUIRE_MESSAGE(reencoding == current, "Error for conversion with value " << value);		
		value += 1;
	} while (value != 0);
}

BOOST_AUTO_TEST_CASE(TestOppositeProd)
{
    StandardPIF<16> minus_one(
			0,
			28,
			1,
			0,
			0
		);
 	for (int16_t i = 0 ; i >= 0 ; ++i) {
        StandardPositEncoding<16> enc{i};
        StandardPositEncoding<16> opposite{i * -1};
		auto value = posit_decoder(enc);
		auto neg_prod = posit_mul(value, minus_one);
		auto op_value = posit_decoder(opposite);
		auto convert_op = PositIF_to_PositProd(op_value);
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
    StandardPositEncoding<16> val{0};
	auto decoded = posit_decoder(val);
	BOOST_REQUIRE_MESSAGE(decoded.getExp() == 0, 
			"Decoded biased exp of 0 should be zero after decoding"
		);
}

BOOST_AUTO_TEST_CASE(TestQuireConvertBack)
{
    StandardQuire<16> quire{0};

	for(uint32_t value = 0; value < (1<<16); value++) {
        StandardPositEncoding<16> valueEncoding{value};
        StandardPIF<16> decoded{valueEncoding};
        StandardPositProd<16> prod{decoded};
        StandardQuire<16> quireConvert = add_sub_quire(quire, prod, 0);
		auto back_convert = quire_to_posit(quireConvert);
		if (decoded.getIsNaR() == 0 and back_convert != decoded) {
			cerr << "=== Original : ===" << endl;
			decoded.printContent();
			cerr << "=== Quire : ===" << endl;
			cerr << quireConvert << endl;
			cerr << "=== Decoded : ===" << endl;
			back_convert.printContent();
		}
		if (decoded.getIsNaR() == 1) {
			BOOST_REQUIRE_MESSAGE(back_convert.getIsNaR() == 1 ,
				"Nar value decoding should be NaR"
			);
		} else {
			BOOST_REQUIRE_MESSAGE(back_convert == decoded,
					"Error for posit with encoding " << value
				);
		}
	}
}

BOOST_AUTO_TEST_CASE(TestSegmentedQuireConvertBack)
{
    StandardSegmentedQuire<16, 64> quire{0};

	for(uint32_t value = 0; value < (1<<16); value++) {
        StandardPositEncoding<16> valueEncoding{value};
        StandardPIF<16> decoded{valueEncoding};
        StandardPositProd<16> prod{valueEncoding};
        StandardSegmentedQuire<16, 64> segmentedQuireConvert = segmented_add_sub_quire(quire, prod, 0);
        StandardQuire<16>  quireConvert = propagateCarries(segmentedQuireConvert);
		auto back_convert = quire_to_posit(quireConvert);
		if (decoded.getIsNaR() == 0 and back_convert != decoded) {
			cerr << "=== Original : ===" << endl;
			decoded.printContent();
			cerr << "=== Segmented Quire : ===" << endl;
			segmentedQuireConvert.printContent();
			cerr << "=== Quire : ===" << endl;
			quireConvert.printContent();
			cerr << "=== Decoded : ===" << endl;
			back_convert.printContent();
		}
		if (decoded.getIsNaR() == 1) {
			BOOST_REQUIRE_MESSAGE(back_convert.getIsNaR() == 1 ,
				"Nar value decoding should be NaR"
			);
		} else {
			BOOST_REQUIRE_MESSAGE(back_convert == decoded,
					"Error for posit with encoding " << value
				);
		}
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

    array<StandardPIF<16>, 4> values = {
        StandardPIF<16>::getMinPos(),
        StandardPIF<16>::getMaxPos(),
        StandardPIF<16>::getMinNeg(),
        StandardPIF<16>::getMaxNeg()
	};

	array<string, 4> names = {
		"MinPos", 
		"MaxPos",
		"MinNeg",
		"MaxNeg"
	};

	for (size_t i = 0 ; i < 4 ; ++i) {
		uint16_t enc = encodings[i];
        StandardPIF<16> & valStat = values[i];
		string& name = names[i];
        auto decoded = posit_decoder(StandardPositEncoding<16>{enc});
		BOOST_REQUIRE_MESSAGE(decoded == valStat, 
				"Error in " << name << " generated value should be identical to decoded one");
	}

}

BOOST_AUTO_TEST_CASE(PositToValueTestPosit16) 
{
	uint16_t i = 0;
	do {
		//cout << "Testing " << i << endl;
        StandardPositEncoding<16> encoding(i);
		auto decoded = posit_decoder(encoding);
		ap_uint<14> complete_sig = decoded.getSignBit().concat(decoded.getSignificand());
		int value = complete_sig.to_int();
		if (decoded.getSignBit() == 1) {
            value = twoCompClean(value, StandardPositDim<16>::WF + 2);
		}
		
        double floatingPointVal = value;

        int exp = decoded.getExp() - StandardPositDim<16>::EXP_BIAS;
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
        StandardPositEncoding<16> encoding(i);
		auto decoded = posit_decoder(encoding);
		auto encoded = posit_encoder(decoded);

		BOOST_REQUIRE_MESSAGE(encoding == encoded,
				"Error : " << i << " is not encoded correctly");
		i += 1;
	} while (i != 0); 
}
