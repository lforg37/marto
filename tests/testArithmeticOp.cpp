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


template<unsigned int N, bool is_signed>
using Wrapper = hint::VivadoWrapper<N, is_signed>;

bool is_mul_ok_p16(uint16_t testval0, uint16_t testval1)
{
	StandardPositEncoding<16, Wrapper> marto_p_0{{testval0}}, marto_p_1{{testval1}};
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
	// ap_uint<16> value10 = 0b1000000000000001;
	// ap_uint<16> value20 = 0b1000000000000010;
	// PositValue<16> a = posit_decoder(value10);
	// PositValue<16> b = posit_decoder(value20);
	// // PositValue<16> b = posit_decoder((ap_uint<16>)0b0110111100000111);

	// PositProd<16> res = posit_mul(a,b);

	// printf("===== a =====\n");
	// a.printContent();
	// printf("===== b =====\n");
	// b.printContent();
	// printf("===== res =====\n");
	// res.printContent();
	// auto resPosit = PositProd_to_PositValue(res);
	// printf("===== res as posit =====\n");
	// resPosit.printContent();

	// ap_uint<16> encoded = posit_encoder(resPosit);
	// printApUint(encoded);

	// posit16_t positValue1 = castP16(value10);
	// posit16_t positValue2 = castP16(value20);
	// posit16_t positMul = p16_mul(positValue1, positValue2);
	// ap_uint<16> softpositMul = (ap_uint<16>) castUI(positMul);
	// printApUint(softpositMul);

	// printf("===== encoding of soft result =====\n");
	// ap_uint<16> value30 = softpositMul;
	// PositValue<16> t = posit_decoder(value30);
	// t.printContent();
	// exit(0);

	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1}<<32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = StandardPositEncoding<16, Wrapper>{{value2}};
		auto decoded2 = StandardPIF<16, Wrapper, true>{value2Encoding};

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, Wrapper>{{value1}};
			auto prod = value1Encoding * value2Encoding;
			//cerr << to_string(static_cast<Wrapper<StandardPositDim<16>::ProdSize, false> >(prod)) << endl;
			StandardPositEncoding<16, Wrapper> encoded{prod};
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positMul = p16_mul(positValue1, positValue2);
			Wrapper<16, false> softpositMul{castUI(positMul)};
			if(!(encoded == softpositMul).template isSet<0>()){
				/*fprintf(stderr, "\n\n\n\n");
				fprintf(stderr, "=== Inputs === \n");
				printApUint(value1Encoding);
				printApUint(value2Encoding);
				fprintf(stderr, "=== Expected result === \n");
				printApUint(softpositMul);
				fprintf(stderr, "=== Computed result === \n");
				printApUint(encoded);
				fprintf(stderr, "Tests Passed: %lu\n", counter);
				*/
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

BOOST_AUTO_TEST_CASE(TestAllSubQuirePosit16, *utf::disabled() * utf::label("long"))
{
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1} << 32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = StandardPositEncoding<16, Wrapper> {{value2}};
		auto decoded2 = posit_decoder(value2Encoding);
		auto prod2 = PositIF_to_PositProd(decoded2);
		auto base_quire = add_sub_quire(StandardQuire<16, Wrapper>{}, prod2, {0});

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, Wrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto prod1 = PositIF_to_PositProd(decoded1);
			auto sub = add_sub_quire(base_quire, prod1, {1});
			//cerr << to_string(static_cast<Wrapper<StandardQuireDim<16>::Size, false> >(sub)) <<
			//	endl;
			auto subval = quire_to_posit(sub);
			//cerr <<
			//	to_string(static_cast<Wrapper<StandardPositDim<16>::ValSize, false> >(subval)) <<
			//	endl;
			auto encoded = posit_encoder(subval);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_sub(positValue2, positValue1);
			Wrapper<16, false> softpositSum{castUI(positSum)};
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
	StandardPositEncoding<16, Wrapper> marto_p_0{{testval0}}, marto_p_1{{testval1}};
	auto decoded0 = posit_decoder(marto_p_0);
	auto decoded1 = posit_decoder(marto_p_1);
	posit16_t positValue0 = castP16(testval0);
	posit16_t positValue1 = castP16(testval1);
	posit16_t positSum = p16_add(positValue0, positValue1);
	auto sum = posit_add(decoded0, decoded1);
	auto final = posit_encoder(sum);
	return final.unravel() == castUI(positSum);
}

BOOST_AUTO_TEST_CASE(TestSomeSumP16)
{
	BOOST_REQUIRE(is_sum_ok_p16(0, 1));
}

BOOST_AUTO_TEST_CASE(TestAllSumPosit16, *utf::disabled() * utf::label("long"))
{
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1}<<32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = StandardPositEncoding<16, Wrapper>{{value2}};
		auto decoded2 = posit_decoder(value2Encoding);

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, Wrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto sum = posit_add(decoded1, decoded2);
			auto encoded = posit_encoder(sum);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_add(positValue1, positValue2);
			Wrapper<16, false> softpositSum {castUI(positSum)};
			if(!(encoded == softpositSum).template isSet<0>()){
				/*fprintf(stderr, "\n\n\n\n");
				fprintf(stderr, "=== Inputs === \n");
				printApUint(value1Encoding);
				printApUint(value2Encoding);
				fprintf(stderr, "=== Expected result === \n");
				printApUint(softpositSum);
				fprintf(stderr, "=== Computed result === \n");
				printApUint(encoded);
				sum.printContent();
				fprintf(stderr, "Tests Passed: %lu\n", counter);
				*/
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
		auto value2Encoding =  - StandardPositEncoding<16, Wrapper>{{value2}};
		// value2Encoding = ~value2Encoding+1;
		auto decoded2 = posit_decoder(value2Encoding);

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, Wrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto sub = posit_add(decoded1, decoded2);
			auto encoded = posit_encoder(sub);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSub = p16_sub(positValue1, positValue2);
			Wrapper<16, false> softpositSum{castUI(positSub)};
			if(!(encoded == softpositSum).template isSet<0>()){
				BOOST_REQUIRE_MESSAGE(false, "Sub of " << value1 << " and " << value2 << " returned " << hint::to_string(encoded) << " while it should have returned " << to_string(softpositSum));
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
		auto value2Encoding = StandardPositEncoding<16, Wrapper>{{value2}};
		auto decoded2 = posit_decoder(value2Encoding);
		auto prod2 = PositIF_to_PositProd(decoded2);
		auto base_quire = segmented_add_sub_quire(StandardSegmentedQuire<16, 16, Wrapper>{}, prod2, {0});
		auto quire = add_sub_quire(StandardQuire<16, Wrapper>{}, prod2, {0});

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, Wrapper>{{value1}};
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
			Wrapper<16, false> softpositSum {{castUI(positSum)}};
			if(!(encoded == softpositSum).isSet<0>()){
				BOOST_REQUIRE_MESSAGE(false, "Sum of " << value1 << " and " << value2 << " returned " << hint::to_string(encoded) << " while it should have returned " << hint::to_string(softpositSum));
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
	auto value2Encoding = StandardPositEncoding<16, Wrapper>{{val2}};
	auto decoded2 = posit_decoder(value2Encoding);
	auto value1Encoding = StandardPositEncoding<16, Wrapper>{{val1}};
	auto decoded1 = posit_decoder(value1Encoding);
	PositIntermediateFormat<16, 1, Wrapper, true> sum = posit_add_in_place(decoded1, decoded2);
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
		auto value2Encoding = StandardPositEncoding<16, Wrapper>{{value2}};
		auto decoded2 = posit_decoder(value2Encoding);

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, Wrapper>{{value1}};
			auto decoded1 = posit_decoder(value1Encoding);
			auto sum = posit_add_in_place(decoded1, decoded2);
			auto encoded = posit_encoder(sum);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_add(positValue1, positValue2);
			Wrapper<16, false> softpositSum {castUI(positSum)};
			// int32_t hard = encoded.unravel();
			// int32_t soft = softpositSum.unravel();
			// int32_t diff = ((hard-soft)<0) ? soft-hard : hard-soft;

			if(!(encoded == softpositSum).template isSet<0>()){
			// if(diff > 1){
				/*fprintf(stderr, "\n\n\n\n");
				fprintf(stderr, "=== Inputs === \n");
				printApUint(value1Encoding);
				printApUint(value2Encoding);
				fprintf(stderr, "=== Expected result === \n");
				printApUint(softpositSum);
				fprintf(stderr, "=== Computed result === \n");
				printApUint(encoded);
				sum.printContent();
				fprintf(stderr, "Tests Passed: %lu\n", counter);
				*/
				BOOST_REQUIRE_MESSAGE(false, "Sum of " << value1 << " and " << value2 << " returned " << to_string(encoded) << " while it should have returned " << to_string(softpositSum));
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

#if 0
BOOST_AUTO_TEST_CASE(TestAllMulPosit16_in_place_rounding, *utf::disabled() * utf::label("long"))
{
	// auto value1Encoding_single = StandardPositEncoding<16, Wrapper>{{1}};
	// auto value2Encoding_single = StandardPositEncoding<16, Wrapper>{{1}};
	// auto decoded1_single = posit_decoder(value1Encoding_single);
	// auto decoded2_single = posit_decoder(value2Encoding_single);
	// auto prod = posit_mul(decoded1_single, decoded2_single);
	// auto mul_in_place = PositProd_to_PositIF_in_place_rounding(prod);
	// auto prod_single = posit_mul(decoded1_single, decoded2_single);
	// auto mul_single = PositProd_to_PositIF(prod_single);
	// auto encoded_in_place_single = posit_encoder(mul_in_place);
	// auto encoded_single = posit_encoder(mul_single);
	// auto decoded_mul_single = posit_decoder(encoded_single);

	// cerr << "value1Encoding_single: " << to_string(value1Encoding_single) << endl;
	// cerr << "value2Encoding_single: " << to_string(value2Encoding_single) << endl;
	// cerr << "decoded1_single: " << to_string(decoded1_single) << endl;
	// cerr << "decoded2_single: " << to_string(decoded2_single) << endl;
	// cerr << "mul_in_place:        " << to_string(mul_in_place) << endl;
	// cerr << "mul_single:          " << to_string(mul_single) << endl;
	// cerr << "decoded_mul_single:  " << to_string(decoded_mul_single) << endl;
	// cerr << "encoded_in_place_single: " << to_string(encoded_in_place_single) << endl;
	// cerr << "encoded_single:          " << to_string(encoded_single) << endl;

	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = uint64_t{1}<<32;
	unsigned int error_counter = 0;
	#pragma omp parallel for
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = StandardPositEncoding<16, Wrapper>{{value2}};
		auto decoded2 = StandardPIF<16, Wrapper, true>{value2Encoding};

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, Wrapper>{{value1}};
			auto decoded1 = StandardPIF<16, Wrapper, true>{value1Encoding};
			auto prod = posit_mul(decoded1, decoded2);
			auto pif = PositProd_to_PositIF_in_place_rounding(prod);
			//cerr << to_string(static_cast<Wrapper<StandardPositDim<16>::ProdSize, false> >(prod)) << endl;
			StandardPositEncoding<16, Wrapper> encoded{pif};
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positMul = p16_mul(positValue1, positValue2);
			Wrapper<16, false> softpositMul{castUI(positMul)};
			if(!(encoded == softpositMul).template isSet<0>()){
				/*fprintf(stderr, "\n\n\n\n");
				fprintf(stderr, "=== Inputs === \n");
				printApUint(value1Encoding);
				printApUint(value2Encoding);
				fprintf(stderr, "=== Expected result === \n");
				printApUint(softpositMul);
				fprintf(stderr, "=== Computed result === \n");
				printApUint(encoded);
				fprintf(stderr, "Tests Passed: %lu\n", counter);
				*/
				BOOST_REQUIRE_MESSAGE(false, "Mul of " << value1 << " and " << value2 << " returned " << to_string(encoded) << " while it should have returned " << to_string(softpositMul));
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
#endif
