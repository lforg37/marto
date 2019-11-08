#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PositToValueTest

#define VIVADO_BACKEND
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>

#include <boost/test/unit_test.hpp>
#ifdef SOFTPOSIT
#include "softposit.h"
#endif

#include <primitives/lzoc_shifter.hpp>
#include <tools/printing.hpp>

using hint::to_string;

#include "posit/posit_dim.hpp"
#include "posit/add_sub_quire.hpp"
#include "posit/quire_to_posit.hpp"

using namespace std;
template<unsigned int N, bool is_signed>
using Wrapper = hint::VivadoWrapper<N, is_signed>;

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
	StandardQuire<16, Wrapper> quire{};
	constexpr unsigned int PIF_SIZE = StandardPositDim<16>::ValSize;
	constexpr unsigned int PROD_SIZE = StandardPositDim<16>::ProdSize;
	constexpr unsigned int QUIRE_SIZE = StandardQuireDim<16>::Size;
	auto minpos = StandardPIF<16, Wrapper>::getMinPos();
	auto prod = posit_mul(minpos, minpos);
	auto quire_conv = add_sub_quire(quire, prod, {0});
	auto decoded = quire_to_posit(quire_conv);
	bool ok = (decoded == minpos).isSet<0>();
	if (!ok) {
		cerr << "minpos:\t" << hint::to_string(static_cast<Wrapper<PIF_SIZE, false> >(minpos)) << endl;
		cerr << "prod:\t" << hint::to_string(static_cast<Wrapper<PROD_SIZE, false> >(prod)) << endl;
		cerr << "quire:\t" << hint::to_string(static_cast<Wrapper<QUIRE_SIZE, false> >(quire_conv)) << endl;
		cerr << "decoded:\t" << hint::to_string(static_cast<Wrapper<PIF_SIZE, false> >(decoded)) << endl;

	}
	BOOST_REQUIRE_MESSAGE(ok,
			"Positive underflow does not returns minpos"
		);

	//Positive sticky
	prod = PositIF_to_PositProd(minpos);
	quire_conv = add_sub_quire(quire_conv, prod, {0});
	decoded = quire_to_posit(quire_conv);
	auto encoding = posit_encoder(decoded);
	auto minpos_enc = posit_encoder(minpos);
	auto res = StandardPIF<16, Wrapper>{{1<<19}};
	BOOST_REQUIRE_MESSAGE((encoding == minpos_enc).isSet<0>(),
			"Min pos * (1+minpos) should return minpos"
		);
	BOOST_REQUIRE_MESSAGE(decoded.getStickyBit().isSet<0>(), "minpos * 1+minpos should set sticky bit");

	//Positive overflow
	auto maxpos = StandardPIF<16, Wrapper>::getMaxPos();
	prod = posit_mul(maxpos, maxpos);
	quire_conv = add_sub_quire(quire, prod, {0});
	decoded = quire_to_posit(quire_conv);
	BOOST_REQUIRE_MESSAGE((decoded == maxpos).isSet<0>(),
			"Positive overflow doesn't return maxpos"
			);

	//Negative underflow
	auto minneg = StandardPIF<16, Wrapper>::getMinPos();
	prod = posit_mul(minpos, minneg);
	quire_conv = add_sub_quire(quire, prod, {0});
	decoded = quire_to_posit(quire_conv);
	BOOST_REQUIRE_MESSAGE((decoded == minneg).isSet<0>(),
			"Negative underflow is not mapped to minneg"
			);

	//Negative overflow
	auto maxneg = StandardPIF<16, Wrapper>::getMaxPos();
	prod = posit_mul(maxpos, maxneg);
	quire_conv = add_sub_quire(quire, prod, {0});
	decoded = quire_to_posit(quire_conv);
	BOOST_REQUIRE_MESSAGE((decoded == maxneg).isSet<0>(),
			"Negative overflow is not mapped to maxneg"
			);
}

BOOST_AUTO_TEST_CASE(PositValueToProd)
{
	uint16_t value = 0;
	constexpr unsigned int PROD_SIZE = StandardPositDim<16>::ProdSize;
	StandardPIF<16, Wrapper> one({0}, {StandardPositDim<16>::EXP_BIAS}, {0}, {1}, {0});
	do {
		PositEncoding<16, 1, Wrapper> posit_encoding{{value}};
		auto posit_val = static_cast<PositIntermediateFormat<16, 1, Wrapper> >(posit_encoding);

		auto posit_prod_direct = static_cast<PositProd<16, 1, Wrapper> >(posit_val);
		auto posit_prod_by_one = posit_mul(posit_val, one);

		if (posit_val.getIsNaR().template isSet<0>()) {
			BOOST_REQUIRE_MESSAGE(posit_prod_by_one.getIsNaR().template isSet<0>() , "Prod by one should be NAR");
			BOOST_REQUIRE_MESSAGE(posit_prod_direct.getIsNaR().template isSet<0>(), "Direct conversion should be NAR");
		} else {
			bool ok = (posit_prod_by_one == posit_prod_direct).isSet<0>();
			if (!ok) {
				cerr << "posit_prod_by_one : " << hint::to_string(static_cast<Wrapper<PROD_SIZE, false> >(posit_prod_by_one)) << endl;
				cerr << "posit_prod_direct : " << hint::to_string(static_cast<Wrapper<PROD_SIZE, false> >(posit_prod_direct)) << endl;
			}
			BOOST_REQUIRE_MESSAGE(ok, "Error for conversion with value " << value);
		}
		value += 1;
	} while (value != 0);
}

BOOST_AUTO_TEST_CASE(PositValueToProdToValue)
{
	uint16_t value = 0;
	do {
		StandardPositEncoding<16, Wrapper> current{{value}};
		auto decoded = static_cast<StandardPIF<16, Wrapper> >(current);
		auto prod = static_cast<StandardPositProd<16, Wrapper> >(current);
		auto casted_val = PositProd_to_PositIF(prod);
		auto reencoding = static_cast<StandardPositEncoding<16, Wrapper> >(casted_val);

		BOOST_REQUIRE_MESSAGE((reencoding == current).isSet<0>(), "Error for conversion with value " << value);
		value += 1;
	} while (value != 0);
}

BOOST_AUTO_TEST_CASE(TestOppositeProd)
{
	StandardPIF<16, Wrapper> minus_one(
			{0},
			{28},
			{1},
			{0},
			{0}
		);
	for (int16_t i = 0 ; i >= 0 ; ++i) {
		StandardPositEncoding<16, Wrapper> enc{{i}};
		StandardPositEncoding<16, Wrapper> opposite{{i * -1}};
		auto value = posit_decoder(enc);
		auto neg_prod = posit_mul(value, minus_one);
		auto op_value = posit_decoder(opposite);
		auto convert_op = PositIF_to_PositProd(op_value);
		BOOST_REQUIRE_MESSAGE(
				(convert_op == neg_prod).isSet<0>(),
				"Error with encoding " << i
			);
	}
}


BOOST_AUTO_TEST_CASE(TestZeroExpZero)
{
	StandardPositEncoding<16, Wrapper> val{{0}};
	Wrapper<StandardPositDim<16>::WE, false> zero{0};
	auto decoded = posit_decoder(val);
	BOOST_REQUIRE_MESSAGE((decoded.getExp() == zero).isSet<0>(),
			"Decoded biased exp of 0 should be zero after decoding"
		);
}

BOOST_AUTO_TEST_CASE(TestQuireConvertBack)
{
	StandardQuire<16, Wrapper> quire{};

	//TODO activvae whole loop
	for(uint32_t value = 0; value < (1<<16); value++) {
		StandardPositEncoding<16, Wrapper> valueEncoding{{value}};
		StandardPIF<16, Wrapper> decoded{valueEncoding};
		StandardPositProd<16, Wrapper> prod{decoded};
		StandardQuire<16, Wrapper> quireConvert = add_sub_quire(quire, prod, {0});
		auto back_convert = quire_to_posit(quireConvert);
		if (decoded.getIsNaR().template isSet<0>()) {
			BOOST_REQUIRE_MESSAGE(back_convert.getIsNaR().template isSet<0>(),
				"Nar value decoding should be NaR"
			);
		} else {
			BOOST_REQUIRE_MESSAGE((back_convert == decoded).isSet<0>(),
					"Error for posit with encoding " << value
				);
		}
	}
}

BOOST_AUTO_TEST_CASE(TestSegmentedQuireConvertBack)
{
	StandardSegmentedQuire<16, 64, Wrapper> quire{};

	//TODO restore loop
	for(uint32_t value = 0; value < (1<<16); value++) {
		StandardPositEncoding<16, Wrapper> valueEncoding{{value}};
		StandardPIF<16, Wrapper> decoded{valueEncoding};
		StandardPositProd<16, Wrapper> prod{valueEncoding};
		StandardSegmentedQuire<16, 64, Wrapper> segmentedQuireConvert = segmented_add_sub_quire(quire, prod, {0});
		//cerr << to_string(static_cast<Wrapper<segmentedQuireConvert.Size, false> >(segmentedQuireConvert)) << endl;
		StandardQuire<16, Wrapper>  quireConvert = propagateCarries(segmentedQuireConvert);
		auto back_convert = quire_to_posit(quireConvert);
		if (decoded.getIsNaR().template isSet<0>()) {
			BOOST_REQUIRE_MESSAGE(back_convert.getIsNaR().template isSet<0>() ,
				"Nar value decoding should be NaR"
			);
		} else {
			BOOST_REQUIRE_MESSAGE((back_convert == decoded).isSet<0>(),
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

	array<StandardPIF<16, Wrapper>, 4> values = {
		StandardPIF<16, Wrapper>::getMinPos(),
		StandardPIF<16, Wrapper>::getMaxPos(),
		StandardPIF<16, Wrapper>::getMinNeg(),
		StandardPIF<16, Wrapper>::getMaxNeg()
	};

	array<string, 4> names = {
		"MinPos",
		"MaxPos",
		"MinNeg",
		"MaxNeg"
	};

	for (size_t i = 0 ; i < 4 ; ++i) {
		uint16_t enc = encodings[i];
		StandardPIF<16, Wrapper> & valStat = values[i];
		string& name = names[i];
		auto decoded = posit_decoder(StandardPositEncoding<16, Wrapper>{{enc}});
		BOOST_REQUIRE_MESSAGE((decoded == valStat).isSet<0>(),
				"Error in " << name << " generated value should be identical to decoded one");
	}

}
/*
BOOST_AUTO_TEST_CASE(PositToValueTestPosit16)
{
	uint16_t i = 0;
	do {
		//cout << "Testing " << i << endl;
		StandardPositEncoding<16, Wrapper> encoding{{i}};
		auto decoded = posit_decoder(encoding);
		auto complete_sig = decoded.getSignBit().concatenate(decoded.getSignificand());
		int value = complete_sig.to_int();
		if (decoded.getSignBit() == 1) {
			value = twoCompClean(value, StandardPositDim<16>::WF + 2);
		}

		double floatingPointVal = value;

		int exp = decoded.getExp() - StandardPositDim<16>::EXP_BIAS;
		double factor = pow(2.0, (exp-12));
		floatingPointVal *= factor;
		//cout << "Decoded value : " << value << endl;

		posit16_t posit_val = castP16(i);
		double soft_posit_val = convertP16ToDouble(posit_val);

		if (decoded.getIsNaR() == 1) {
			BOOST_REQUIRE_MESSAGE(isnan(soft_posit_val),
					"Error in conversion : decoding of value " << i <<
					"gives NaN when softposit value gives " << soft_posit_val
					);
		} else {
			//cout << "Soft Posit decoded value : " << soft_posit_val << endl << endl;
			BOOST_REQUIRE_MESSAGE(
					soft_posit_val == floatingPointVal,
					"Error in conversion : decoding of value " << i <<
					" gives " << floatingPointVal << " instead of " << soft_posit_val );
		}
		i += 1;
	}while(i != 0);
}*/
#endif

BOOST_AUTO_TEST_CASE(PositEncoder)
{
	// Should be tested on 32 bits as well when possible
	uint16_t i = 0;
	do {
		StandardPositEncoding<16, Wrapper> encoding{{i}};
		auto decoded = posit_decoder(encoding);
		auto encoded = posit_encoder(decoded);

		BOOST_REQUIRE_MESSAGE((encoding == encoded).isSet<0>(),
				"Error : " << i << " is not encoded correctly");
		i += 1;
	} while (i != 0);
}
