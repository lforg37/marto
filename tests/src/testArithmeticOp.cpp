#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PositArithmeticOps

#ifndef VIVADO_BACKEND
#define VIVADO_BACKEND
#endif

#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>

#include <boost/test/unit_test.hpp>
#include "softposit.h"

#include "posit/posit_dim.hpp"
#include "kulisch/kulisch_dim.hpp"
#include "posit/add_sub_quire.hpp"
//#include "tools/.hpp"
#include "ieeefloats/ieee_adder.hpp"

#include "hint.hpp"
#include "tools/printing.hpp"

#include <omp.h>

using namespace std;
using hint::to_string;
namespace utf = boost::unit_test;


using hint::VivadoWrapper;

bool is_mul_ok_p16(uint16_t testval0, uint16_t testval1)
{
	StandardPositEncoding<16, VivadoWrapper> marto_p_0{{testval0}}, marto_p_1{{testval1}};
	auto decoded0 = posit_decoder(marto_p_0);
	auto decoded1 = posit_decoder(marto_p_1);
	posit16_t positValue0 = castP16(testval0);
	posit16_t positValue1 = castP16(testval1);
	posit16_t positProd = p16_mul(positValue0, positValue1);
	auto prod = posit_mul(decoded0, decoded1);
	auto final = posit_encoder(PositProd_to_PositIF(prod));
	bool res = final.unravel() == castUI(positProd);
	if (!res) {
		std::cerr << "Expecting " << castUI(positProd) << std::endl;
	}
	return res;
}

BOOST_AUTO_TEST_CASE(TestSomeMulPosit16) {
	BOOST_REQUIRE(is_mul_ok_p16(1, 4096));
	BOOST_REQUIRE(is_mul_ok_p16(2, 24576));
	BOOST_REQUIRE(is_mul_ok_p16(1, 8192));
	BOOST_REQUIRE(is_mul_ok_p16(1, 49152));
	BOOST_REQUIRE(is_mul_ok_p16(5, 57344));
	BOOST_REQUIRE(is_mul_ok_p16(5, 40960));
	BOOST_REQUIRE(is_mul_ok_p16(49152, 47104));
	BOOST_REQUIRE(is_mul_ok_p16(32769, 57344));
}

BOOST_AUTO_TEST_CASE(TestAllMulPosit16, *utf::disabled() * utf::label("long"))
{
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1}<<32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = StandardPositEncoding<16, VivadoWrapper>{{value2}};
		auto decoded2 = StandardPIF<16, VivadoWrapper, true>{value2Encoding};

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, VivadoWrapper>{{value1}};
			auto prod = value1Encoding * value2Encoding;
			StandardPositEncoding<16, VivadoWrapper> encoded{prod};
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positMul = p16_mul(positValue1, positValue2);
			VivadoWrapper<16, false> softpositMul{castUI(positMul)};
			if(!(encoded == softpositMul).template isSet<0>()){
				BOOST_REQUIRE_MESSAGE(false, "Mul of " << value1 << " and " << value2 << " returned a problematic value while it should have returned " << to_string(softpositMul));
			}
		}
		if(((value2%100) == 0) and (value2 != 0)){
			#pragma omp atomic
			counter+=(100*(1<<16));
			#pragma omp critical
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(TOTAL_TESTS)*100, counter,TOTAL_TESTS);
		}
		error_counter = 0;
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)\n", static_cast<double>(TOTAL_TESTS)/static_cast<double>(TOTAL_TESTS)*100, TOTAL_TESTS,TOTAL_TESTS);
}

BOOST_AUTO_TEST_CASE(TestAllMulPosit8)
{
	for(uint32_t value2 = 0; value2 < (1<<8); value2++){
		auto value2Encoding = StandardPositEncoding<8, VivadoWrapper>{{value2}};
		auto decoded2 = StandardPIF<8, VivadoWrapper, true>{value2Encoding};

		for(uint32_t value1 = 0; value1 < (1<<8); value1++){
			auto value1Encoding = StandardPositEncoding<8, VivadoWrapper>{{value1}};
			auto prod = value1Encoding * value2Encoding;
			StandardPositEncoding<8, VivadoWrapper> encoded{prod};
			posit8_t positValue1 = castP8(value1);
			posit8_t positValue2 = castP8(value2);
			posit8_t positMul = p8_mul(positValue1, positValue2);
			VivadoWrapper<8, false> softpositMul{castUI(positMul)};
			BOOST_REQUIRE_MESSAGE((encoded == softpositMul).unravel(), "Mul of " << value1 << " and " << value2 << " returned a problematic value while it should have returned " << to_string(softpositMul));
		}
	}
}

BOOST_AUTO_TEST_CASE(TestAllSubQuirePosit16, *utf::disabled() * utf::label("long"))
{
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1} << 32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = StandardPositEncoding<16, VivadoWrapper> {{value2}};
		auto decoded2 = posit_decoder(value2Encoding);
		auto prod2 = PositIF_to_PositProd(decoded2);
		auto base_quire = add_sub_quire(StandardQuire<16, VivadoWrapper>{}, prod2, {0});

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, VivadoWrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto prod1 = PositIF_to_PositProd(decoded1);
			auto sub = add_sub_quire(base_quire, prod1, {1});
			//cerr << to_string(static_cast<VivadoWrapper<StandardQuireDim<16>::Size, false> >(sub)) <<
			//	endl;
			auto subval = quire_to_posit(sub);
			//cerr <<
			//	to_string(static_cast<VivadoWrapper<StandardPositDim<16>::ValSize, false> >(subval)) <<
			//	endl;
			auto encoded = posit_encoder(subval);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_sub(positValue2, positValue1);
			VivadoWrapper<16, false> softpositSum{castUI(positSum)};
			auto ok = (encoded == softpositSum);
			bool res = ok.isSet<0>();
			if(not res){
				BOOST_REQUIRE_MESSAGE(false, "Sum of " << value1 << " and " << value2 << " returned a wrong value while it should have returned an other value");
			}
		}
		if(((value2%20) == 0) and (value2 != 0)){
			#pragma omp atomic
			counter+=(20*(1<<16));
			#pragma omp critical
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%%  (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(TOTAL_TESTS)*100, counter,TOTAL_TESTS);
		}
		error_counter = 0;
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%%  (%lu\t/%lu)\n", static_cast<double>(TOTAL_TESTS)/static_cast<double>(TOTAL_TESTS)*100, TOTAL_TESTS,TOTAL_TESTS);
}

bool is_sum_ok_p16(uint16_t testval0, uint16_t testval1)
{
	StandardPositEncoding<16, VivadoWrapper> marto_p_0{{testval0}}, marto_p_1{{testval1}};
	auto decoded0 = posit_decoder(marto_p_0);
	auto decoded1 = posit_decoder(marto_p_1);
	posit16_t positValue0 = castP16(testval0);
	posit16_t positValue1 = castP16(testval1);
	posit16_t positSum = p16_add(positValue0, positValue1);
	auto sum = posit_add(decoded0, decoded1);
	auto final = posit_encoder(sum);
	return final.unravel() == castUI(positSum);
}

bool is_sum_ok_p8(uint8_t testval0, uint8_t testval1)
{
	StandardPositEncoding<8, VivadoWrapper> marto_p_0{{testval0}}, marto_p_1{{testval1}};
	auto decoded0 = posit_decoder(marto_p_0);
	auto decoded1 = posit_decoder(marto_p_1);
	posit8_t positValue0 = castP8(testval0);
	posit8_t positValue1 = castP8(testval1);
	posit8_t positSum = p8_add(positValue0, positValue1);
	auto sum = posit_add(decoded0, decoded1);
	auto final = posit_encoder(sum);
	return final.unravel() == castUI(positSum);
}

BOOST_AUTO_TEST_CASE(TestSomeSumP16)
{
	BOOST_REQUIRE(is_sum_ok_p16(0, 1));
	BOOST_REQUIRE(is_sum_ok_p16((1 << 13) - 1, ((1 << 13) - 1) * -1));
	BOOST_REQUIRE(is_sum_ok_p16(-1*((1 << 13) - 1), (1 << 13) - 1));
}

BOOST_AUTO_TEST_CASE(TestSomeSumP8)
{
	BOOST_REQUIRE(is_sum_ok_p8(127, 127));
	BOOST_REQUIRE(is_sum_ok_p8(129, 129));
}

BOOST_AUTO_TEST_CASE(TestAllSumPosit8)
{
	for(uint32_t value2 = 0; value2 < (1<<8); value2++){
		auto value2Encoding = StandardPositEncoding<8, VivadoWrapper>{{value2}};
		auto decoded2 = posit_decoder(value2Encoding);
		for(uint32_t value1 = 0; value1 < (1<<8); value1++){
			auto value1Encoding = StandardPositEncoding<8, VivadoWrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto sum = posit_add(decoded1, decoded2);
			auto encoded = posit_encoder(sum);
			posit8_t positValue1 = castP8(value1);
			posit8_t positValue2 = castP8(value2);
			posit8_t positSum = p8_add(positValue1, positValue2);
			VivadoWrapper<8, false> softpositSum {castUI(positSum)};
			BOOST_REQUIRE_MESSAGE((encoded == softpositSum).template isSet<0>(),
								  "Sum of " << value1 << " and " << value2 << " returned a wrong value while it should have returned " << to_string(softpositSum));
		}
	}
}

BOOST_AUTO_TEST_CASE(TestAllSumPosit16, *utf::disabled() * utf::label("long"))
{
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1}<<32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = StandardPositEncoding<16, VivadoWrapper>{{value2}};
		auto decoded2 = posit_decoder(value2Encoding);

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, VivadoWrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto sum = posit_add(decoded1, decoded2);
			auto encoded = posit_encoder(sum);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_add(positValue1, positValue2);
			VivadoWrapper<16, false> softpositSum {castUI(positSum)};
			if(!(encoded == softpositSum).template isSet<0>()){
				BOOST_REQUIRE_MESSAGE(false, "Sum of " << value1 << " and " << value2 << " returned a wrong value while it should have returned " << to_string(softpositSum));
			}
		}
		if(((value2%100) == 0) and (value2 != 0)){
			#pragma omp atomic
			counter+=(100*(1<<16));
			#pragma omp critical
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%%  (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(TOTAL_TESTS)*100, counter,TOTAL_TESTS);
		}
		error_counter = 0;
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%%  (%lu\t/%lu)\n", static_cast<double>(TOTAL_TESTS)/static_cast<double>(TOTAL_TESTS)*100, TOTAL_TESTS,TOTAL_TESTS);
}


BOOST_AUTO_TEST_CASE(TestAllSubPosit16, *utf::disabled() * utf::label("long"))
{
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1} << 32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding =  - StandardPositEncoding<16, VivadoWrapper>{{value2}};
		// value2Encoding = ~value2Encoding+1;
		auto decoded2 = posit_decoder(value2Encoding);

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, VivadoWrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto sub = posit_add(decoded1, decoded2);
			auto encoded = posit_encoder(sub);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSub = p16_sub(positValue1, positValue2);
			VivadoWrapper<16, false> softpositSum{castUI(positSub)};
			if(!(encoded == softpositSum).template isSet<0>()){
				BOOST_REQUIRE_MESSAGE(false, "Sub of " << value1 << " and " << value2 << " returned " << hint::to_string(static_cast<VivadoWrapper<16, false> const &>(encoded)) << " while it should have returned " << to_string(softpositSum));
			}
		}
		if(((value2%100) == 0) and (value2 != 0)){
			#pragma omp atomic
			counter+=(100*(1<<16));
			#pragma omp critical
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%%  (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(TOTAL_TESTS)*100, counter,TOTAL_TESTS);
		}
		error_counter = 0;
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%%  (%lu\t/%lu)\n", static_cast<double>(TOTAL_TESTS)/static_cast<double>(TOTAL_TESTS)*100, TOTAL_TESTS,TOTAL_TESTS);
}

BOOST_AUTO_TEST_CASE(TestStaticDivide)
{
	BOOST_REQUIRE_MESSAGE((Static_Ceil_Div<4,2>::val == 2), "Error with value 4/2.");
	BOOST_REQUIRE_MESSAGE((Static_Ceil_Div<5,2>::val == 3), "Error with value 5/2.");
}

BOOST_AUTO_TEST_CASE(TestAllSegmentedSubQuirePosit16, *utf::disabled() * utf::label("long"))
{
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1}<<32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = StandardPositEncoding<16, VivadoWrapper>{{value2}};
		auto decoded2 = posit_decoder(value2Encoding);
		auto prod2 = PositIF_to_PositProd(decoded2);
		auto base_quire = segmented_add_sub_quire(StandardSegmentedQuire<16, 16, VivadoWrapper>{}, prod2, {0});
		auto quire = add_sub_quire(StandardQuire<16, VivadoWrapper>{}, prod2, {0});

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, VivadoWrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto prod1 = PositIF_to_PositProd(decoded1);
			auto sub = segmented_add_sub_quire(base_quire, prod1, {1});
			//cerr << to_string(sub) << endl;
			auto sub_quire = add_sub_quire(quire, prod1, {1});
			auto propagation = propagateCarries(sub);
			auto subval = quire_to_posit(propagation);
			auto encoded = posit_encoder(subval);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_sub(positValue2, positValue1);
			VivadoWrapper<16, false> softpositSum {{castUI(positSum)}};
			if(!(encoded == softpositSum).isSet<0>()){
				BOOST_REQUIRE_MESSAGE(false, "Sum of " << value1 << " and " << value2 << " returned " << hint::to_string(static_cast<VivadoWrapper<16, false> const &>(encoded)) << " while it should have returned " << hint::to_string(softpositSum));
			}
		}
		if(((value2%20) == 0) and (value2 != 0)){
			#pragma omp atomic
			counter+=(20*(1<<16));
			#pragma omp critical
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%%  (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(TOTAL_TESTS)*100, counter,TOTAL_TESTS);
		}
		error_counter = 0;
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)\n", static_cast<double>(TOTAL_TESTS)/static_cast<double>(TOTAL_TESTS)*100, TOTAL_TESTS,TOTAL_TESTS);
}

inline bool test_some_inplace_sum_p16(uint16_t val1, uint16_t val2)
{
	auto value2Encoding = StandardPositEncoding<16, VivadoWrapper>{{val2}};
	auto decoded2 = posit_decoder(value2Encoding);
	auto value1Encoding = StandardPositEncoding<16, VivadoWrapper>{{val1}};
	auto decoded1 = posit_decoder(value1Encoding);
	PositIntermediateFormat<16, 1, VivadoWrapper, true> sum = posit_add_in_place(decoded1, decoded2);
	auto encoded = posit_encoder(sum);
	posit16_t positValue1 = castP16(val1);
	posit16_t positValue2 = castP16(val2);
	posit16_t positSum = p16_add(positValue1, positValue2);
	uint16_t positSumInt = castUI(positSum);
	bool res = (positSumInt == encoded.unravel());
	if (not res) {
		auto sumunrounded = posit_add(decoded1, decoded2);
		auto directencode = posit_encoder(sumunrounded);
		std::cerr << (directencode.unravel() == positSumInt) << std::endl;
	}
	return res;
}

BOOST_AUTO_TEST_CASE(TestSomeSumPosit16_in_place_rounding)
{
	BOOST_REQUIRE(test_some_inplace_sum_p16(16383, 16384));
}

BOOST_AUTO_TEST_CASE(TestAllSumPosit16_in_place_rounding, *utf::disabled() * utf::label("long"))
{
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1}<<32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = StandardPositEncoding<16, VivadoWrapper>{{value2}};
		auto decoded2 = posit_decoder(value2Encoding);

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, VivadoWrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto sum = posit_add_in_place(decoded1, decoded2);
			auto encoded = posit_encoder(sum);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_add(positValue1, positValue2);
			VivadoWrapper<16, false> softpositSum {castUI(positSum)};

			if(!(encoded == softpositSum).template isSet<0>()){
			BOOST_REQUIRE_MESSAGE(false, "Sum of " << value1 << " and " << value2 << " returned " << to_string(static_cast<VivadoWrapper<16, false> const &>(encoded)) << " while it should have returned " << to_string(softpositSum));
			}
		}
		if(((value2%100) == 0) and (value2 != 0)){
			#pragma omp atomic
			counter+=(100*(1<<16));
			#pragma omp critical
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%%  (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(TOTAL_TESTS)*100, counter,TOTAL_TESTS);
		}
		error_counter = 0;
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%%  (%lu\t/%lu)\n", static_cast<double>(TOTAL_TESTS)/static_cast<double>(TOTAL_TESTS)*100, TOTAL_TESTS,TOTAL_TESTS);
}

BOOST_AUTO_TEST_CASE(TestAllSumPosit8_in_place_rounding)
{
	for(uint32_t value2 = 0; value2 < (1<<8); value2++){
		auto value2Encoding = StandardPositEncoding<8, VivadoWrapper>{{value2}};
		auto decoded2 = posit_decoder(value2Encoding);

		for(uint32_t value1 = 0; value1 < (1<<8); value1++){
			auto value1Encoding = StandardPositEncoding<8, VivadoWrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto sum = posit_add_in_place(decoded1, decoded2);
			auto encoded = posit_encoder(sum);
			posit8_t positValue1 = castP8(value1);
			posit8_t positValue2 = castP8(value2);
			posit8_t positSum = p8_add(positValue1, positValue2);
			VivadoWrapper<8, false> softpositSum {castUI(positSum)};
			BOOST_REQUIRE_MESSAGE((encoded == softpositSum).unravel(), "Sum of " << value1 << " and " << value2 <<
								  " returned " << to_string(static_cast<VivadoWrapper<8, false> const &>(encoded)) << " while it should have returned " << to_string(softpositSum));
		}
	}
}

BOOST_AUTO_TEST_CASE(TestAllMulPosit16_in_place_rounding, *utf::disabled() * utf::label("long"))
{
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1}<<32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = StandardPositEncoding<16, VivadoWrapper>{{value2}};
		auto decoded2 = StandardPIF<16, VivadoWrapper, true>{value2Encoding};

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, VivadoWrapper>{{value1}};
			auto decoded1 = StandardPIF<16, VivadoWrapper, true>{value1Encoding};
			auto prod = posit_mul(decoded1, decoded2);
			auto pif = PositProd_to_PositIF_in_place_rounding(prod);
			//cerr << to_string(static_cast<VivadoWrapper<StandardPositDim<16>::ProdSize, false> >(prod)) << endl;
			StandardPositEncoding<16, VivadoWrapper> encoded{pif};
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positMul = p16_mul(positValue1, positValue2);
			VivadoWrapper<16, false> softpositMul{castUI(positMul)};
			if(!(encoded == softpositMul).template isSet<0>()){
				BOOST_REQUIRE_MESSAGE(false, "Mul of " << value1 << " and " << value2 << " returned " << to_string(static_cast<VivadoWrapper<16, false> const &>(encoded)) << " while it should have returned " << to_string(softpositMul));
			}
		}
		if(((value2%100) == 0) and (value2 != 0)){
			#pragma omp atomic
			counter+=(100*(1<<16));
			#pragma omp critical
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(TOTAL_TESTS)*100, counter,TOTAL_TESTS);
		}
		error_counter = 0;
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)\n", static_cast<double>(TOTAL_TESTS)/static_cast<double>(TOTAL_TESTS)*100, TOTAL_TESTS,TOTAL_TESTS);
}


BOOST_AUTO_TEST_CASE(TestAllMulPosit8_in_place_rounding)
{
	for(uint32_t value2 = 0; value2 < (1<<8); value2++){
		auto value2Encoding = StandardPositEncoding<8, VivadoWrapper>{{value2}};
		auto decoded2 = StandardPIF<8, VivadoWrapper, true>{value2Encoding};

		for(uint32_t value1 = 0; value1 < (1<<8); value1++){
			auto value1Encoding = StandardPositEncoding<8, VivadoWrapper>{{value1}};
			auto decoded1 = StandardPIF<8, VivadoWrapper, true>{value1Encoding};
			auto prod = posit_mul(decoded1, decoded2);
			auto pif = PositProd_to_PositIF_in_place_rounding(prod);
			StandardPositEncoding<8, VivadoWrapper> encoded{pif};
			posit8_t positValue1 = castP8(value1);
			posit8_t positValue2 = castP8(value2);
			posit8_t positMul = p8_mul(positValue1, positValue2);
			VivadoWrapper<8, false> softpositMul{castUI(positMul)};
			BOOST_REQUIRE_MESSAGE((encoded == softpositMul).unravel(), "Mul of " << value1 << " and " << value2 << " returned " <<
								  to_string(static_cast<VivadoWrapper<8, false> const &>(encoded)) << " while it should have returned " << to_string(softpositMul));
		}
	}
}
