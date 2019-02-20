#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE PositArithmeticOps

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/test/unit_test.hpp>
#include "softposit.h"

#include "add_sub_quire.hpp"
#include "lzoc_shifter.hpp"
#include "posit_add.hpp"
#include "posit_decoder.hpp"
#include "posit_dim.hpp"
#include "posit_encoder.hpp"
#include "posit_mul.hpp"
#include "quire_to_posit.hpp"
#include "value_prod_conversions.hpp"

#include <omp.h>

using namespace std;
namespace utf = boost::unit_test;

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
	uint64_t TOTAL_TESTS = (((uint64_t)1)<<32);
	unsigned int error_counter = 0;
	#pragma omp parallel for 
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = PositEncoding<16> (value2);
		auto decoded2 = posit_decoder(value2Encoding);

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = PositEncoding<16> (value1);
			auto decoded1 = posit_decoder(value1Encoding);
			auto sum = posit_mul(decoded1, decoded2);
			auto sumAsPositValue = PositProd_to_PositValue(sum);
			auto encoded = posit_encoder(sumAsPositValue);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positMul = p16_mul(positValue1, positValue2);
			ap_uint<16> softpositMul = (ap_uint<16>) castUI(positMul);
			if(!(encoded == softpositMul)){
				fprintf(stderr, "\n\n\n\n");
				fprintf(stderr, "=== Inputs === \n");
				printApUint(value1Encoding);
				printApUint(value2Encoding);
				fprintf(stderr, "=== Expected result === \n");
				printApUint(softpositMul);
				fprintf(stderr, "=== Computed result === \n");
				printApUint(encoded);
				fprintf(stderr, "Tests Passed: %lu\n", counter);

				BOOST_REQUIRE_MESSAGE(false, "Mul of " << value1 << " and " << value2 << " returned " << (unsigned int)encoded << " while it should have returned " << (unsigned int)softpositMul);
			}
		}
		if(((value2%100) == 0) and (value2 != 0)){
			#pragma omp atomic
			counter+=(100*(1<<16));
			#pragma omp critical
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f\% (%lu\t/%lu)", ((double)counter/(double)TOTAL_TESTS)*100, counter,TOTAL_TESTS);
		}
		error_counter = 0;
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f\% (%lu\t/%lu)\n", ((double)TOTAL_TESTS/(double)TOTAL_TESTS)*100, TOTAL_TESTS,TOTAL_TESTS);
}

BOOST_AUTO_TEST_CASE(TestAllSubQuirePosit16, *utf::disabled() * utf::label("long")) 
{
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = (((uint64_t)1)<<32);
	unsigned int error_counter = 0;
	#pragma omp parallel for 
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = PositEncoding<16> (value2);
		auto decoded2 = posit_decoder(value2Encoding);
		auto prod2 = PositValue_to_PositProd(decoded2);
		auto base_quire = add_sub_quire(Quire<16>{0}, prod2, 0);

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = PositEncoding<16> (value1);
			auto decoded1 = posit_decoder(value1Encoding);
			auto prod1 = PositValue_to_PositProd(decoded1);
			auto sub = add_sub_quire(base_quire, prod1, 1);
			auto subval = quire_to_posit(sub);
			auto encoded = posit_encoder(subval);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_sub(positValue2, positValue1);
			ap_uint<16> softpositSum = (ap_uint<16>) castUI(positSum);
			if(!(encoded == softpositSum)){
				fprintf(stderr, "\n\n\n\n");
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

				BOOST_REQUIRE_MESSAGE(false, "Sum of " << value1 << " and " << value2 << " returned " << (unsigned int)encoded << " while it should have returned " << (unsigned int)softpositSum);
			}
		}
		if(((value2%100) == 0) and (value2 != 0)){
			#pragma omp atomic
			counter+=(100*(1<<16));
			#pragma omp critical
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f\%  (%lu\t/%lu)", ((double)counter/(double)TOTAL_TESTS)*100, counter,TOTAL_TESTS);
		}
		error_counter = 0;
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f\%  (%lu\t/%lu)\n", ((double)TOTAL_TESTS/(double)TOTAL_TESTS)*100, TOTAL_TESTS,TOTAL_TESTS);
}

BOOST_AUTO_TEST_CASE(TestAllSumPosit16, *utf::disabled() * utf::label("long")) 
{

	// ap_uint<16> value10 = 0b1111111111111111;
	// ap_uint<16> value20 = 0b0000000000000001;
	// PositValue<16> a = posit_decoder(value10);
	// PositValue<16> b = posit_decoder(value20);
	// // PositValue<16> b = posit_decoder((ap_uint<16>)0b0110111100000111);

	// PositValue<16> res = posit_add(a,b);  
	
	// printf("===== a =====\n");
	// a.printContent();
	// printf("===== b =====\n");
	// b.printContent();
	// printf("===== res =====\n");
	// res.printContent();

	// ap_uint<16> encoded = posit_encoder(res);
	// printApUint(encoded);

	// posit16_t positValue1 = castP16(value10);
	// posit16_t positValue2 = castP16(value20);
	// posit16_t positSum = p16_add(positValue1, positValue2);
	// ap_uint<16> softpositSum = (ap_uint<16>) castUI(positSum);
	// printApUint(softpositSum);

	// printf("===== encoding of soft result =====\n");
	// ap_uint<16> value30 = softpositSum;
	// PositValue<16> t = posit_decoder(value30);
	// t.printContent();
	// exit(0);

	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = (((uint64_t)1)<<32);
	unsigned int error_counter = 0;
	#pragma omp parallel for 
	for(uint32_t value2 = 0; value2 < (1<<16); value2++){
		auto value2Encoding = PositEncoding<16> (value2);
		auto decoded2 = posit_decoder(value2Encoding);

		for(uint32_t value1 = 0; value1 < (1<<16); value1++){
			auto value1Encoding = PositEncoding<16> (value1);
			auto decoded1 = posit_decoder(value1Encoding);
			auto sum = posit_add(decoded1, decoded2);
			auto encoded = posit_encoder(sum);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_add(positValue1, positValue2);
			ap_uint<16> softpositSum = (ap_uint<16>) castUI(positSum);
			if(!(encoded == softpositSum)){
				fprintf(stderr, "\n\n\n\n");
				fprintf(stderr, "=== Inputs === \n");
				printApUint(value1Encoding);
				printApUint(value2Encoding);
				fprintf(stderr, "=== Expected result === \n");
				printApUint(softpositSum);
				fprintf(stderr, "=== Computed result === \n");
				printApUint(encoded);
				sum.printContent();
				fprintf(stderr, "Tests Passed: %d\n", counter);

				BOOST_REQUIRE_MESSAGE(false, "Sum of " << value1 << " and " << value2 << " returned " << (unsigned int)encoded << " while it should have returned " << (unsigned int)softpositSum);
			}
		}
		if(((value2%100) == 0) and (value2 != 0)){
			#pragma omp atomic
			counter+=(100*(1<<16));
			#pragma omp critical
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f\%  (%lu\t/%lu)", ((double)counter/(double)TOTAL_TESTS)*100, counter,TOTAL_TESTS);
		}
		error_counter = 0;
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f\%  (%lu\t/%lu)\n", ((double)TOTAL_TESTS/(double)TOTAL_TESTS)*100, TOTAL_TESTS,TOTAL_TESTS);
}

BOOST_AUTO_TEST_CASE(TestAddSubQuire)
{
	for(uint32_t value = 0; value < (1<<16); value++) {
		auto valueEncoding = PositEncoding<16> (value);
		auto decoded = posit_decoder(valueEncoding);
	}
}
