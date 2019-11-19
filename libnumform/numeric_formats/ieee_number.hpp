#ifndef IEEENUMBER_HPP
#define IEEENUMBER_HPP

#include <gmpxx.h>
#include <mpfr.h>


namespace libnumform {
	/**
	 * An abstraction of IEEE FP numbers.
	 * User can convert between value (mpfr) and representation (mpz_class)
	 */
	template<size_t WE, size_t WF>
	class IEEENumber
	{
	public:

		/** Several possible special values */
		typedef enum {
			plusInfty,     /**< A positive infinity with random non-zero exponent and fraction bits  */
			minusInfty,    /**< A negative infinity with random non-zero exponent and fraction bits  */
			plusZero,      /**< A positive zero */
			minusZero,     /**< A negative zero  */
			NaN,            /**< Not A Number */
			smallestSubNormal,      /**< The smallest subnormal IEEENumber*/
			greatestSubNormal,      /**< The greatest subnormal IEEENumber*/
			minusGreatestSubNormal,      /**< The negative greatest subnormal IEEENumber*/
			smallestNormal,                 /**< The smallest normal IEEENumber*/
			greatestNormal                  /**< The greatest normal IEEENumber*/
		} SpecialValue;

		/**
		 * Constructs a new IEEENumber.
		 */
		IEEENumber();

		/**
		 * Constructs a new IEEENumber.
		 * @param v a special value
		 */
		IEEENumber(SpecialValue v);

		/**
		 * Constructs a new initialised IEEENumber.
		 * @param m the initial value.
		 * @param ternaryRoundingInfo used to avoid double rounding.
		 */
		IEEENumber(mpfr_t m,  int ternaryRoundInfo=0);

		/**
		 * Constructs a new initialised IEEENumber.
		 * @param z the initial value, given as an mpz holding the bits of the Number.
		 */
		IEEENumber(mpz_class z);

		/**
		 * Retrieves the significant.
		 * @return Returns an mpz_class, representing the
		 * VHDL signal of the mantissa, without leading 1.
		 */
		mpz_class getMantissaSignalValue();


		/**
		 * Converts the currently stored IEEENumber to an mpfr_t
		 * @return a reference to a mpfr containing the represented value
		 */
		 void getMPFR(mpfr_t mp_);


		/**
		 * imports a MPFR value into an IEEENumber
		 * @param[in] m the value
		 * @param  ternaryRoundInfo a ternary rounding value rememberin the rounding history, see mpfr documenation
		 */
		//void setMPFR(mpfr_t const m, int ternaryRoundInfo);

		void setMPFR(mpfr_t m, int ternaryRoundInfo);


		/**
		 * Assignes a signal value. Converts the signal value to the
		 * relevant IEEENumber fields.
		 * @param s the signal value to assign.
		 */
		IEEENumber &operator=(mpz_class s);

		/**
		 * Retrieved the binary signal representation of this floating point.
		 * @return a mpz_class.
		 */
		mpz_class getSignalValue();

		/**
		 * Returns wE and wF.
		 * @param[out] wE_ exponent precision
		 * @param[out] wF_ fraction precision
		 */
		void getPrecision(int &wE, int &wF);

		class IEEEContextManager
		{
			public:
				IEEEContextManager()
				{
					_oldmin = mpfr_get_emin();
					_oldmax = mpfr_get_emax();

					constexpr int _expshiftmax = -1 * (1 << (WE - 1));
					constexpr int _expshiftmaxwfoffest = _expshiftmax - static_cast<int>(WF);
					constexpr int emin = _expshiftmaxwfoffest + 3; // -1024 - 52 + 3
					mpfr_set_emin (emin);
					// The formula should evaluatempfr_t mp to 1024 for doubles, see MPFR doc;
					constexpr int emax = (1<<(WE-1));
					mpfr_set_emax (emax);
				}

				~IEEEContextManager()
				{
					mpfr_set_emax(_oldmax);
					mpfr_set_emin(_oldmin);
				}

				IEEEContextManager(IEEEContextManager const &) = delete;
				IEEEContextManager(IEEEContextManager &&) = delete;

				IEEEContextManager& operator=(IEEEContextManager const &) = delete;
				IEEEContextManager& operator=(IEEEContextManager &&) = delete;

			private:
				mpfr_exp_t _oldmin;
				mpfr_exp_t _oldmax;
		};


	private:
		/** The value of the sign field */
		mpz_class sign;

		/** The value of the exponent field */
		mpz_class exponent;

		/** The value of the mantissa field  */
		mpz_class mantissa;

		static constexpr int _iWE = static_cast<int>(WE);
		static constexpr int _iWF = static_cast<int>(WF);
	};

}

#include "ieee_number.tpp"
#endif

