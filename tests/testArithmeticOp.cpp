#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE PositArithmeticOps

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/test/unit_test.hpp>
#include "softposit.h"

#include "posit_decoder.hpp"
#include "posit_encoder.hpp"
#include "posit_dim.hpp"
#include "posit_add.hpp"
#include "lzoc_shifter.hpp"

#include <omp.h>

using namespace std;

BOOST_AUTO_TEST_CASE(TestAllSumPosit8) 
{
	// ap_uint<16> value10 = 0b1111111111111111;
	// ap_uint<16> value20 = 0b0000000000000000;
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

	// uint16_t value1 = 0;
	uint16_t value2 = 0;
	uint64_t counter = 0;
	uint64_t TOTAL_TESTS = (((uint64_t)1)<<32);
	unsigned int error_counter = 0;
	#pragma omp parallel for private(value2)
	for(uint16_t value1 = 0; value1 <((1<<16)-1); value1++){
		auto value2Encoding = PositEncoding<16> (value2);
		auto decoded2 = posit_decoder(value2Encoding);

		do {
			auto value1Encoding = PositEncoding<16> (value1);
			auto decoded1 = posit_decoder(value1Encoding);
			auto sum = posit_add(decoded1, decoded2);
			auto encoded = posit_encoder(sum);
			posit16_t positValue1 = castP16(value1);
			posit16_t positValue2 = castP16(value2);
			posit16_t positSum = p16_add(positValue1, positValue2);
			ap_uint<16> softpositSum = (ap_uint<16>) castUI(positSum);
			// ap_uint<16> abs_diff;
			// if(encoded>softpositSum){
			// 	abs_diff = encoded-softpositSum;
			// }
			// else{
			// 	abs_diff = softpositSum-encoded;
			// }
			if(!(encoded == softpositSum)){
			// if(log2((double)abs_diff)>1){
				// error_counter++;
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
			counter++;
			value1++;
		} while (value1 != 0);
		value2++;
		if((value2%100) == 0){
			#pragma omp atomic
			counter+=100;
			#pragma omp critical
			fprintf(stderr, "steps completed : %d\t/%lu\n", counter, TOTAL_TESTS);
		}
		error_counter = 0;
	}while (value2 != 0);
}
