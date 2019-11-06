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
namespace utf = boost::unit_test;


template<unsigned int N, bool is_signed>
using Wrapper = hint::VivadoWrapper<N, is_signed>;

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
		auto decoded2 = StandardPIF<16, Wrapper>{value2Encoding};

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = StandardPositEncoding<16, Wrapper>{{value1}};
			StandardPositEncoding<16, Wrapper> encoded{value1Encoding * value2Encoding};
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positMul = p16_mul(positValue1, positValue2);
			ap_uint<16> softpositMul{castUI(positMul)};
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
				BOOST_REQUIRE_MESSAGE(false, "Mul of " << value1 << " and " << value2 << " returned a problematic value while it should have returned " << static_cast<unsigned int>(softpositMul));
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
			auto subval = quire_to_posit(sub);
			auto encoded = posit_encoder(subval);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_sub(positValue2, positValue1);
			Wrapper<16, false> softpositSum{castUI(positSum)};
			auto ok = (encoded == softpositSum);
			bool res = ok.isSet<0>();
			if(not res){
				/*fprintf(stderr, "\n\n\n\n");
				fprintf(stderr, "=== Inputs === \n");
				cerr << "  ";
				printApUint(value2Encoding);
				cerr << "- ";
				printApUint(value1Encoding);
				cerr << "=== Quire ===" << endl;
				printApUint(sub);
				fprintf(stderr, "=== Expected result === \n");
				printApUint(softpositSum);
				fprintf(stderr, "=== Computed result === \n");
				printApUint(encoded);
				subval.printContent();
				fprintf(stderr, "Tests Passed: %lu\n", counter);
*/
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
			ap_uint<16> softpositSum {castUI(positSum)};
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
				BOOST_REQUIRE_MESSAGE(false, "Sum of " << value1 << " and " << value2 << " returned a wrong value while it should have returned " << (unsigned int)softpositSum);
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
			ap_uint<16> softpositSum{castUI(positSub)};
			if(!(encoded == softpositSum).template isSet<0>()){
				BOOST_REQUIRE_MESSAGE(false, "Sub of " << value1 << " and " << value2 << " returned " << hint::to_string(encoded) << " while it should have returned " << static_cast<unsigned int>(softpositSum));
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

/*
BOOST_AUTO_TEST_CASE(TestIEEAddWE3WF4Vivado)
{
	constexpr unsigned int WF = 4;
	constexpr unsigned int WE = 3;

	uint64_t nbTests = 1 << ((1+WF+WE) * 2);

	auto filepath = TESTDATA_ROOT"/ieeeadder/test_WE3_WF4.input";
	ifstream testfile;
	testfile.open(filepath);
	bitset<1+WE+WF> i0;
	bitset<1+WE+WF> i1;
	string tmp;
	bitset<1+WE+WF> res;
	BOOST_REQUIRE_MESSAGE(!testfile.fail(), "The test file " << filepath << " cannot be opened.");
	for (uint64_t curtest = 0; curtest < nbTests ; ++curtest) {
		testfile >> i0 >> i1 >> tmp >> res;

		uint64_t i0num = static_cast<uint64_t>(i0.to_ulong());
		uint64_t i1num = static_cast<uint64_t>(i1.to_ulong());
		uint64_t resnum = static_cast<uint64_t>(res.to_ulong());

		IEEENumber<WE, WF, hint::VivadoWrapper> i0ieee{{i0num}};
		IEEENumber<WE, WF, hint::VivadoWrapper> i1ieee{{i1num}};

		auto result = ieee_add_sub_impl(i0ieee, i1ieee);

		uint64_t convertedRes = result.to_uint();

		//cerr << "IN0 : " << endl << i0gmp.get_str(2) << endl;
		//cerr << "IN1 : " << endl << i1gmp.get_str(2) << endl;
		//cerr << "EXPECTED : " << endl << mpzres.get_str(2) << endl;
		//cerr << "GOT : " << endl << result.unravel().get_str(2) << endl;

		BOOST_REQUIRE_MESSAGE(resnum == convertedRes, "Error for iteration " << curtest);
	}
}
*/
