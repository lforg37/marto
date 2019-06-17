#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE KulischAcc

#include <boost/test/unit_test.hpp>
#define AP_INT_MAX_W 4500


#include "marto.hpp"
#include <stdlib.h>
#include <mpfr.h>

namespace utf = boost::unit_test;

#define N 32
#define SEGMENT 32

typedef union {
  unsigned int i;
  float f;
} float_to_fix;

BOOST_AUTO_TEST_CASE(TestRand2CK1_32, *utf::disabled() * utf::label("long")) 
{
	int SIZE = 10000000;
	float *tab1, *tab2, sign1, sign2;
	ap_uint<N> *tab1_ap, *tab2_ap;

	tab1 = (float*) malloc(sizeof(float)*SIZE);
	tab2 = (float*) malloc(sizeof(float)*SIZE);
	tab1_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);
	tab2_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);

	for(int i=0; i<SIZE; i++){
		tab1[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		tab2[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		sign1 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		sign2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		if(((int) sign1) == 0){
			tab1[i] = -tab1[i]; 
		}
		if(((int) sign2) == 0){
			tab2[i] = -tab2[i]; 
		}

	}

	// tab1[0] = 293.992; 
	// tab2[0] = 0.19882; 


	for(int i=0; i<SIZE; i++){
		float_to_fix conv;
		conv.f = tab1[i] ;
		tab1_ap[i] = conv.i;
		conv.f = tab2[i] ;
		tab2_ap[i] = conv.i;
	};


	mpfr_t mpfr_acc;
	mpfr_init2(mpfr_acc,1000);
	mpfr_t current_mpfr_value1, current_mpfr_value2;
	mpfr_init2(current_mpfr_value1,1000);
	mpfr_init2(current_mpfr_value2,1000);
 	mpfr_set_d(mpfr_acc, 0.0, GMP_RNDZ);
		
	KulischAcc<N> acc(0);
	ap_uint<32> res;
	float mpfr_sum;
	for(int i=0; i<SIZE; i++){
	// for(int i=0; i<1; i++){
		// fprintf(stderr, "tab1 : %f, tab2: %f\n", tab1[i], tab2[i]);
		auto prod = exact_prod<N>(tab1_ap[i], tab2_ap[i]);
		acc = add_2CK1<N>(acc, prod);
		
		mpfr_set_flt(current_mpfr_value1, tab1[i], GMP_RNDZ);
		mpfr_set_flt(current_mpfr_value2, tab2[i], GMP_RNDZ);
		mpfr_mul(current_mpfr_value1, current_mpfr_value1, current_mpfr_value2, GMP_RNDZ);
		mpfr_add(mpfr_acc, mpfr_acc, current_mpfr_value1, GMP_RNDZ);
		
		res = acc_IEEE_rounding<N>(acc);
		float_to_fix conv_r;
		conv_r.i = res;
		float res_f = conv_r.f;
		mpfr_sum = mpfr_get_flt(mpfr_acc, MPFR_RNDN);
		
		if(mpfr_sum != res_f){
			fprintf(stderr, "Result : %1.15f\n", res_f);
			printApUint(res);
			fprintf(stderr, "mpfr : %1.15f\n", mpfr_sum);
			conv_r.f = mpfr_sum;
			ap_uint<32> mpfr_bits = conv_r.i;
			printApUint(mpfr_bits);

			free(tab1);
			free(tab2);
			free(tab1_ap);
			free(tab2_ap);
			BOOST_REQUIRE_MESSAGE(false, "Accumulator holding incorrect result");
		}
		// fprintf(stderr, "%d\n", i);
	}	

	free(tab1);
	free(tab2);
	free(tab1_ap);
	free(tab2_ap);
}


BOOST_AUTO_TEST_CASE(TestRandSMK1_32, *utf::disabled() * utf::label("long")) 
{
	int SIZE = 10000000;
	float *tab1, *tab2, sign1, sign2;
	ap_uint<N> *tab1_ap, *tab2_ap;

	tab1 = (float*) malloc(sizeof(float)*SIZE);
	tab2 = (float*) malloc(sizeof(float)*SIZE);
	tab1_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);
	tab2_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);

	for(int i=0; i<SIZE; i++){
		tab1[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		tab2[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		sign1 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		sign2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		if(((int) sign1) == 0){
			tab1[i] = -tab1[i]; 
		}
		if(((int) sign2) == 0){
			tab2[i] = -tab2[i]; 
		}

	}

	// tab1[0] = 293.992; 
	// tab2[0] = 0.19882; 


	for(int i=0; i<SIZE; i++){
		float_to_fix conv;
		conv.f = tab1[i] ;
		tab1_ap[i] = conv.i;
		conv.f = tab2[i] ;
		tab2_ap[i] = conv.i;
	};


	mpfr_t mpfr_acc;
	mpfr_init2(mpfr_acc,1000);
	mpfr_t current_mpfr_value1, current_mpfr_value2;
	mpfr_init2(current_mpfr_value1,1000);
	mpfr_init2(current_mpfr_value2,1000);
 	mpfr_set_d(mpfr_acc, 0.0, GMP_RNDZ);
		




	SignedKulischAcc<N> acc(0);
	ap_uint<32> res;
	float mpfr_sum;
	for(int i=0; i<SIZE; i++){
		// fprintf(stderr, "tab1 : %f, tab2: %f\n", tab1[i], tab2[i]);
		auto prod = exact_prod<N>(tab1_ap[i], tab2_ap[i]);
		acc = add_SMK1<N>(acc, prod);
		
		mpfr_set_flt(current_mpfr_value1, tab1[i], GMP_RNDZ);
		mpfr_set_flt(current_mpfr_value2, tab2[i], GMP_RNDZ);
		mpfr_mul(current_mpfr_value1, current_mpfr_value1, current_mpfr_value2, GMP_RNDZ);
		mpfr_add(mpfr_acc, mpfr_acc, current_mpfr_value1, GMP_RNDZ);
		
		res = signed_acc_IEEE_rounding<N>(acc);
		float_to_fix conv_r;
		conv_r.i = res;
		float res_f = conv_r.f;
		mpfr_sum = mpfr_get_flt(mpfr_acc, MPFR_RNDN);
		
		if(mpfr_sum != res_f){
			fprintf(stderr, "Result : %1.15f\n", res_f);
			printApUint(res);
			fprintf(stderr, "mpfr : %1.15f\n", mpfr_sum);
			conv_r.f = mpfr_sum;
			ap_uint<32> mpfr_bits = conv_r.i;
			printApUint(mpfr_bits);

			free(tab1);
			free(tab2);
			free(tab1_ap);
			free(tab2_ap);
			BOOST_REQUIRE_MESSAGE(false, "Accumulator holding incorrect result");
		}
		// fprintf(stderr, "%d\n", i);
	}	

	free(tab1);
	free(tab2);
	free(tab1_ap);
	free(tab2_ap);
}


BOOST_AUTO_TEST_CASE(TestRand_segemented_2CK1_32, *utf::disabled() * utf::label("long")) 
{
	int SIZE = 10000000;
	float *tab1, *tab2, sign1, sign2;
	ap_uint<N> *tab1_ap, *tab2_ap;

	tab1 = (float*) malloc(sizeof(float)*SIZE);
	tab2 = (float*) malloc(sizeof(float)*SIZE);
	tab1_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);
	tab2_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);

	for(int i=0; i<SIZE; i++){
		tab1[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		tab2[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		sign1 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		sign2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		if(((int) sign1) == 0){
			tab1[i] = -tab1[i]; 
		}
		if(((int) sign2) == 0){
			tab2[i] = -tab2[i]; 
		}

	}

	// tab1[0] = 293.992; 
	// tab2[0] = 0.19882; 


	for(int i=0; i<SIZE; i++){
		float_to_fix conv;
		conv.f = tab1[i] ;
		tab1_ap[i] = conv.i;
		conv.f = tab2[i] ;
		tab2_ap[i] = conv.i;
	};


	mpfr_t mpfr_acc;
	mpfr_init2(mpfr_acc,1000);
	mpfr_t current_mpfr_value1, current_mpfr_value2;
	mpfr_init2(current_mpfr_value1,1000);
	mpfr_init2(current_mpfr_value2,1000);
 	mpfr_set_d(mpfr_acc, 0.0, GMP_RNDZ);

	acc_segmented_2CK1<N, SEGMENT> acc(KulischAcc<N>(0));
	KulischAcc<N> acc2(0); 
	ap_uint<32> res;
	float mpfr_sum;
	for(int i=0; i<SIZE; i++){
		// fprintf(stderr, "tab1 : %f, tab2: %f\n", tab1[i], tab2[i]);
		auto prod = exact_prod<N>(tab1_ap[i], tab2_ap[i]);
		acc = add_segmented_2CK1<N, SEGMENT>(acc, prod);
		acc2 = add_2CK1<N>(acc2, prod);
		// acc.printContent();
		KulischAcc<N> accprop = propagate_carries_segmented_2CK1(acc);
		
		// printApUint(accprop);
		// printApUint(acc2);

		mpfr_set_flt(current_mpfr_value1, tab1[i], GMP_RNDZ);
		mpfr_set_flt(current_mpfr_value2, tab2[i], GMP_RNDZ);
		mpfr_mul(current_mpfr_value1, current_mpfr_value1, current_mpfr_value2, GMP_RNDZ);
		mpfr_add(mpfr_acc, mpfr_acc, current_mpfr_value1, GMP_RNDZ);
		
		res = acc_IEEE_rounding<N>(accprop);
		float_to_fix conv_r;
		conv_r.i = res;
		float res_f = conv_r.f;
		mpfr_sum = mpfr_get_flt(mpfr_acc, MPFR_RNDN);
		
		if(mpfr_sum != res_f){
			fprintf(stderr, "Result : %1.15f\n", res_f);
			printApUint(res);
			fprintf(stderr, "mpfr : %1.15f\n", mpfr_sum);
			conv_r.f = mpfr_sum;
			ap_uint<32> mpfr_bits = conv_r.i;
			printApUint(mpfr_bits);

			free(tab1);
			free(tab2);
			free(tab1_ap);
			free(tab2_ap);
			BOOST_REQUIRE_MESSAGE(false, "Accumulator holding incorrect result");
		}
		// fprintf(stderr, "%d\n", i);
	}	

	free(tab1);
	free(tab2);
	free(tab1_ap);
	free(tab2_ap);
}

BOOST_AUTO_TEST_CASE(TestRand2CK3_32, *utf::disabled() * utf::label("long")) 
{
	int SIZE = 10000000;
	float *tab1, *tab2, sign1, sign2;
	ap_uint<N> *tab1_ap, *tab2_ap;

	tab1 = (float*) malloc(sizeof(float)*SIZE);
	tab2 = (float*) malloc(sizeof(float)*SIZE);
	tab1_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);
	tab2_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);

	for(int i=0; i<SIZE; i++){
		tab1[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		tab2[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		sign1 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		sign2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		if(((int) sign1) == 0){
			tab1[i] = -tab1[i]; 
		}
		if(((int) sign2) == 0){
			tab2[i] = -tab2[i]; 
		}

	}

	// tab1[0] = 293.992; 
	// tab2[0] = 0.19882; 


	for(int i=0; i<SIZE; i++){
		float_to_fix conv;
		conv.f = tab1[i] ;
		tab1_ap[i] = conv.i;
		conv.f = tab2[i] ;
		tab2_ap[i] = conv.i;
	};


	mpfr_t mpfr_acc;
	mpfr_init2(mpfr_acc,1000);
	mpfr_t current_mpfr_value1, current_mpfr_value2;
	mpfr_init2(current_mpfr_value1,1000);
	mpfr_init2(current_mpfr_value2,1000);
 	mpfr_set_d(mpfr_acc, 0.0, GMP_RNDZ);
		
	acc_2CK3<N, SEGMENT> acc(KulischAcc<N>(0));
	KulischAcc<N> acc2(0); 
	ap_uint<32> res;
	float mpfr_sum;
	for(int i=0; i<SIZE; i++){
		// fprintf(stderr, "tab1 : %f, tab2: %f\n", tab1[i], tab2[i]);
		auto prod = exact_prod<N>(tab1_ap[i], tab2_ap[i]);
		acc = add_2CK3<N, SEGMENT>(acc, prod);
		acc2 = add_2CK1<N>(acc2, prod);
		// acc.printContent();
		KulischAcc<N> accprop = propagate_carries_2CK3(acc);
		
		// printApUint(accprop);
		// printApUint(acc2);

		mpfr_set_flt(current_mpfr_value1, tab1[i], GMP_RNDZ);
		mpfr_set_flt(current_mpfr_value2, tab2[i], GMP_RNDZ);
		mpfr_mul(current_mpfr_value1, current_mpfr_value1, current_mpfr_value2, GMP_RNDZ);
		mpfr_add(mpfr_acc, mpfr_acc, current_mpfr_value1, GMP_RNDZ);
		
		res = acc_IEEE_rounding<N>(accprop);
		float_to_fix conv_r;
		conv_r.i = res;
		float res_f = conv_r.f;
		mpfr_sum = mpfr_get_flt(mpfr_acc, MPFR_RNDN);
		
		if(mpfr_sum != res_f){
			fprintf(stderr, "Result : %1.15f\n", res_f);
			printApUint(res);
			fprintf(stderr, "mpfr : %1.15f\n", mpfr_sum);
			conv_r.f = mpfr_sum;
			ap_uint<32> mpfr_bits = conv_r.i;
			printApUint(mpfr_bits);

			free(tab1);
			free(tab2);
			free(tab1_ap);
			free(tab2_ap);
			BOOST_REQUIRE_MESSAGE(false, "Accumulator holding incorrect result");
		}
		// fprintf(stderr, "%d\n", i);
	}	

	free(tab1);
	free(tab2);
	free(tab1_ap);
	free(tab2_ap);
}



BOOST_AUTO_TEST_CASE(TestRandSMK3_32, *utf::disabled() * utf::label("long")) 
{
	int SIZE = 10000000;
	float *tab1, *tab2, sign1, sign2;
	ap_uint<N> *tab1_ap, *tab2_ap;

	tab1 = (float*) malloc(sizeof(float)*SIZE);
	tab2 = (float*) malloc(sizeof(float)*SIZE);
	tab1_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);
	tab2_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);

	for(int i=0; i<SIZE; i++){
		tab1[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		tab2[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		sign1 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		sign2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		if(((int) sign1) == 0){
			tab1[i] = -tab1[i]; 
		}
		if(((int) sign2) == 0){
			tab2[i] = -tab2[i]; 
		}

	}

	// tab1[0] = 293.992; 
	// tab2[0] = 0.19882; 


	for(int i=0; i<SIZE; i++){
		float_to_fix conv;
		conv.f = tab1[i] ;
		tab1_ap[i] = conv.i;
		conv.f = tab2[i] ;
		tab2_ap[i] = conv.i;
	};


	mpfr_t mpfr_acc;
	mpfr_init2(mpfr_acc,1000);
	mpfr_t current_mpfr_value1, current_mpfr_value2;
	mpfr_init2(current_mpfr_value1,1000);
	mpfr_init2(current_mpfr_value2,1000);
 	mpfr_set_d(mpfr_acc, 0.0, GMP_RNDZ);

	acc_SMK3<N, SEGMENT> acc(KulischAcc<N>(0));
	KulischAcc<N> acc2(0); 
	ap_uint<32> res;
	float mpfr_sum;
	for(int i=1; i<SIZE; i++){
		//fprintf(stderr, "tab1 : %f, tab2: %f\n", tab1[i], tab2[i]);
		auto prod = exact_prod<N>(tab1_ap[i], tab2_ap[i]);
		acc = add_SMK3<N, SEGMENT>(acc, prod);
		acc2 = add_2CK1<N>(acc2, prod);
		//acc.printContent();
		KulischAcc<N> accprop = propagate_carries_SMK3(acc);
		
		//printApUint(accprop);
		//printApUint(acc2);

		mpfr_set_flt(current_mpfr_value1, tab1[i], GMP_RNDZ);
		mpfr_set_flt(current_mpfr_value2, tab2[i], GMP_RNDZ);
		mpfr_mul(current_mpfr_value1, current_mpfr_value1, current_mpfr_value2, GMP_RNDZ);
		mpfr_add(mpfr_acc, mpfr_acc, current_mpfr_value1, GMP_RNDZ);
		
		res = acc_IEEE_rounding<N>(accprop);
		float_to_fix conv_r;
		conv_r.i = res;
		float res_f = conv_r.f;
		mpfr_sum = mpfr_get_flt(mpfr_acc, MPFR_RNDN);
		
		if(mpfr_sum != res_f){
			fprintf(stderr, "Result : %1.15f\n", res_f);
			printApUint(res);
			fprintf(stderr, "mpfr : %1.15f\n", mpfr_sum);
			conv_r.f = mpfr_sum;
			ap_uint<32> mpfr_bits = conv_r.i;
			printApUint(mpfr_bits);

			free(tab1);
			free(tab2);
			free(tab1_ap);
			free(tab2_ap);
			BOOST_REQUIRE_MESSAGE(false, "Accumulator holding incorrect result");
		}
		// fprintf(stderr, "%d\n", i);
	}	

	free(tab1);
	free(tab2);
	free(tab1_ap);
	free(tab2_ap);
}
