/*
  IEEE-compatible floating-point numbers for FloPoCo

  Author: F. de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2010.
  Kindly borroWEd from FloPoCo
  All rights reserved.

  */
#include <random>

namespace libnumform {
	template<size_t WE, size_t WF>
	IEEENumber<WE, WF>::IEEENumber()
	{
		static_assert (WE < 30, "Exponents larger than 30 bits are not supported");
	}

	template<size_t WE, size_t WF>
	IEEENumber<WE, WF>::IEEENumber(mpfr_t mp_, int ternaryRoundInfo)
		: IEEENumber{}
	{
		setMPFR(mp_, ternaryRoundInfo);
	}

	template<size_t WE, size_t WF>
	IEEENumber<WE, WF>::IEEENumber(SpecialValue v)
		: IEEENumber{}
	{
		switch(v)  {
		case plusInfty:
			sign = 0;
			exponent = (mpz_class{1} << WE) -1;
			mantissa = 0;
			break;
		case minusInfty:
			sign = 1;
			exponent = (mpz_class{1} << WE) -1;
			mantissa = 0;
			break;
		case plusZero:
			sign = 0;
			exponent = 0;
			mantissa = 0;
			break;
		case minusZero:
			sign = 1;
			exponent = 0;
			mantissa = 0;
			break;
		case NaN:
			sign = 0;
			exponent = (mpz_class{1} << WE) -1;
			mantissa = 1;
			break;
		case smallestSubNormal:
			sign = 0;
			exponent = 0;
			mantissa = mpz_class{1};
			break;
		case greatestSubNormal:
			sign = 0;
			exponent = 0;
			mantissa = (mpz_class{1} << WF) -1;
			break;
		case minusGreatestSubNormal:
			sign = 1;
			exponent = 0;
			mantissa = (mpz_class{1} << WF) -1;
			break;
		case smallestNormal:
			sign = 0;
			exponent = mpz_class{1};
			mantissa = 0;
			break;
		case greatestNormal:
			sign = 0;
			exponent = (mpz_class{1} << WE) -2;
			mantissa = (mpz_class{1} << WF) -1;
			break;
		}
	}

	template<size_t WE, size_t WF>
	IEEENumber<WE, WF>::IEEENumber(mpz_class z)
		: IEEENumber{}
	{
		operator=(z);
	}

	/**
	template<size_t WE, size_t WF>
	void IEEENumber<WE, WF>::setMPFR(mpfr_t const mp_, int ternaryRoundInfo){
		mpfr_t mp; // will hold a rounded copy of the number

		// emin and emax are specified for a mantissa in (0.5, 1)
		// The formula should evaluate to -1073 for doubles, see MPFR doc;

		mpfr_init2(mp, 1+WF);
		mpfr_set(mp, mp_, MPFR_RNDN);
		mpfr_subnormalize (mp, ternaryRoundInfo, MPFR_RNDN);
		if (mpfr_nan_p(mp))	{
			sign = 0;
			exponent = (1<<WE)-1;
			mantissa = (1<<WF)-1; // qNaN
		}
		else {
			// all the other values are signed
			sign = mpfr_signbit(mp) == 0 ? 0 : 1;

			if (mpfr_inf_p(mp))	{
				exponent = (1<<WE)-1;
				mantissa = 0;
			}
			else {

				if (mpfr_zero_p(mp)) {
					exponent = 0;
					mantissa = 0;
				}
				else{
					mpfr_abs(mp, mp, MPFR_RNDN);
					mp_exp_t exp = mpfr_get_exp(mp)-1;

					//cout << "exp=" << exp <<endl;
					if(exp + ((1<<(WE-1))-1) <=0) {			// subnormal
						// TODO manage double rounding to subnormals
						exponent=0;

						mpfr_mul_2si(mp, mp, _iWF-1+((1<<(_iWE-1))-1), MPFR_RNDN);
						mpfr_get_z(mantissa.get_mpz_t(), mp,  MPFR_RNDN);
						//cout << "subnormal! " << WF + (exp + ((1<<(WE-1))-1)) << " mantissa=" << mantissa << endl;
					}
					else { // Normal number
						mpfr_div_2si(mp, mp, exp, MPFR_RNDN); // exact operation
						mpfr_sub_ui(mp, mp, 1, MPFR_RNDN);    // exact operation
						mpfr_mul_2si(mp, mp, _iWF, MPFR_RNDN);  // exact operation
						mpfr_get_z(mantissa.get_mpz_t(), mp,  MPFR_RNDN); // exact operation


						// Due to rounding, the mantissa might overflow (i.e. become bigger
						// then WE expect).
						if (mantissa == mpz_class(1) << WF)
						{
							exp++;
							mantissa = 0;
						}

						if (mantissa >= mpz_class(1) << WF)
							throw std::string("Mantissa is too big after conversion to VHDL signal.");
						if (mantissa < 0)
							throw std::string("Mantissa is negative after conversion to VHDL signal.");

						exponent = exp + ((1<<(_iWE-1))-1);

						if (exponent >= (1<<_iWE))
						{
							exponent = (1<<_iWE) -1;
							mantissa = 0;
						}
					}
				}
			}

		}
		mpfr_clear(mp);
	}*/

	template<size_t WE, size_t WF>
	void IEEENumber<WE, WF>::setMPFR(mpfr_t mp, int ternaryRoundInfo){
		mpfr_subnormalize (mp, ternaryRoundInfo, MPFR_RNDN);
		/* NaN */
		if (mpfr_nan_p(mp))	{
			sign = 0;
			exponent = (1<<WE)-1;
			mantissa = (1<<WF)-1; // qNaN
		}
		else {
			// all the other values are signed
			sign = mpfr_signbit(mp) == 0 ? 0 : 1;

			/* Inf */
			if (mpfr_inf_p(mp))	{
				exponent = (1<<WE)-1;
				mantissa = 0;
			}
			else {

				/* Zero */
				if (mpfr_zero_p(mp)) {
					exponent = 0;
					mantissa = 0;
				}
				else{

					/* Normal and subnormal numbers */
					mpfr_abs(mp, mp, MPFR_RNDN);

					/* Get exponent
							 * mpfr_get_exp() return exponent for significant in [1/2,1)
							 * but WE use [1,2). Hence the -1.
							 */
					mp_exp_t exp = mpfr_get_exp(mp)-1;

					//cout << "exp=" << exp <<endl;
					if(exp + ((1<<(WE-1))-1) <=0) {			// subnormal
						// TODO manage double rounding to subnormals
						exponent=0;
						/* Extract mantissa */
						mpfr_mul_2si(mp, mp, _iWF-1+((1<<(_iWE-1))-1), MPFR_RNDN);
						mpfr_get_z(mantissa.get_mpz_t(), mp,  MPFR_RNDN);
						//cout << "subnormal! " << WF + (exp + ((1<<(WE-1))-1)) << " mantissa=" << mantissa << endl;
					}
					else { // Normal number
						/* Extract mantissa */
						mpfr_div_2si(mp, mp, exp, MPFR_RNDN); // exact operation
						mpfr_sub_ui(mp, mp, 1, MPFR_RNDN);    // exact operation
						mpfr_mul_2si(mp, mp, _iWF, MPFR_RNDN);  // exact operation
						mpfr_get_z(mantissa.get_mpz_t(), mp,  MPFR_RNDN); // exact operation


						// Due to rounding, the mantissa might overflow (i.e. become bigger
						// then WE expect).
						if (mantissa == mpz_class(1) << WF)
						{
							exp++;
							mantissa = 0;
						}

						if (mantissa >= mpz_class(1) << WF)
							throw std::string("Mantissa is too big after conversion to VHDL signal.");
						if (mantissa < 0)
							throw std::string("Mantissa is negative after conversion to VHDL signal.");

						/* Bias  exponent */
						exponent = exp + ((1<<(_iWE-1))-1);

						/* Handle overflow */
						if (exponent >= (1<<_iWE))
						{
							exponent = (1<<_iWE) -1;
							mantissa = 0;
						}
					}
				}
			}
		}
	}


	template<size_t WE, size_t WF>
	mpz_class IEEENumber<WE, WF>::getSignalValue()
	{
		/* Sanity checks */
		if ((sign != 0) && (sign != 1))
			throw std::string("IEEENumber::getSignal: sign is invalid.");
		if ((exponent < 0) || (exponent >= (1<<WE)))
			throw std::string("IEEENumber::getSignal: exponent is invalid.");
		if ((mantissa < 0) || (mantissa >= (mpz_class(1)<<WF)))
			throw std::string("IEEENumber::getSignal: mantissa is invalid.");
		return ((( sign << WE) + exponent) << WF) + mantissa;
	}

	template<size_t WE, size_t WF>
	void IEEENumber<WE, WF>::getMPFR(mpfr_t mp)
	{

		/* NaN */
		if ( (exponent==((1<<WE)-1)) && mantissa!=0 )
		{
			mpfr_set_nan(mp);
			return;
		}

		/* Infinity */
		if ((exponent==((1<<WE)-1)) && mantissa==0)	{
			mpfr_set_inf(mp, (sign == 1) ? -1 : 1);
			return;
		}

		/* Zero and subnormal numbers */
		if (exponent==0)	{
			mpfr_set_z(mp, mantissa.get_mpz_t(), MPFR_RNDN);
			mpfr_div_2si(mp, mp, WF + ((1<<(WE-1))-2), MPFR_RNDN);
			// Sign
			if (sign == 1)
				mpfr_neg(mp, mp, MPFR_RNDN);
			return;
		} // TODO Check it works with signed zeroes

		/* „Normal” numbers
				 * mp = (-1) * (1 + (mantissa / 2^WF)) * 2^unbiased_exp
				 * unbiased_exp = exp - (1<<(WE-1)) + 1
				 */
		mpfr_set_prec(mp, WF+1);
		mpfr_set_z(mp, mantissa.get_mpz_t(), MPFR_RNDN);
		mpfr_div_2si(mp, mp, WF, MPFR_RNDN);
		mpfr_add_ui(mp, mp, 1, MPFR_RNDN);

		mp_exp_t exp = exponent.get_si();
		exp -= ((1<<(WE-1))-1);
		mpfr_mul_2si(mp, mp, exp, MPFR_RNDN);

		// Sign
		if (sign == 1)
			mpfr_neg(mp, mp, MPFR_RNDN);
	}

	template<size_t WE, size_t WF>
	IEEENumber<WE, WF>& IEEENumber<WE, WF>::operator=(mpz_class s)
	{
		//		cerr << "s=" << s << endl;
		mantissa = s & ((mpz_class{1} << WF) - 1);
		s >>= WF;
		exponent = s & ((mpz_class{1} << WE) - 1);
		s >>= WE;
		sign = s & mpz_class{1};

		if (s > 1)
			throw std::string("IEEENumber::operator= s is bigger than expected.");

		return *this;
	}
}
