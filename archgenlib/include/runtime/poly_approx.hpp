#ifndef RUNTIME_POLY_APPROX_HPP
#define RUNTIME_POLY_APPROX_HPP

#include <optional>
#include <string>
#include <iostream>

#include <sollya.h>
#include <gmpxx.h>

#include "fixedpoint/fixedpoint.hpp"
#include "runtime/sollya_fix_function.hpp"

namespace archgenlib{
    /// Vastly inspired from flopoco code
 	/** The BasicPolyApprox object builds and maintains a machine-efficient polynomial approximation to a fixed-point function over some interval
			Fixed point, hence only absolute errors/accuracy targets.

			There are two families of constructors:
			one that inputs target accuracy and computes degree and approximation error (calling the buildApproxFromTargetAccuracy method)
			one that inputs degree and computes approximation error (calling the buildApproxFromDegreeAndLSBs method).

			The first is useful in standalone, or to evaluate degree/MSBs etc.
			The second is needed in the typical case of a domain split, where the degree is determined when determining the split.

			Sketch of the algorithm for  buildApproxFromTargetAccuracy:
				guessDegree gives a tentative degree.
				target_accuracy defines the best-case LSB of the constant part of the polynomial.
				if  addGuardBitsToConstant, we add g=ceil(log2(degree+1)) bits to the LSB of the constant:
				this provides a bit of freedom to fpminimax, for free in terms of evaluation.
				And we call fpminimax with that, and check what it provided.
				Since x is in [0,1], the MSB of coefficient i is then exactly the MSB of a_i.x^i
				For the same reason, all the coefficients should have the same target LSB

			No domain splitting or other range reduction here: these should be managed in different classes
			No good handling of functions with zero coefficients for now.
			- a function with zero constant should be transformed into a "proper" one outside this class.
			   Example: sin, log(1+x)
		  - a function with odd or even Taylor should be transformed as per the Muller book.

			To implement a generic approximator we will need to lift these restrictions, but it is unclear that it needs to be in this class.
			A short term TODO is to detect such cases.
	*/

	class BasicPolyApprox {
	public:
		/** A minimal constructor that inputs a sollya_obj_t function, a degree and the weight of the LSBs.
				This one is mostly for "internal" use by classes that compute the degree separately, e.g. PiecewisePolyApprox
				@param fS: defines the function, as a sollya_obj_t
				@param degree: degree of the polynomial
				@param LSB: weight of the coefficients
		 */
		BasicPolyApprox(SollyaFunction fS, int degree, int LSBOut);

		/** the degree of this polynomial approximation */ 
		int getDegree();

		/** retrieve the bound on approximation error for this polynomial approximation */ 
		double getApproxErrorBound();

		/** retrieve the i-th coefficient  */ 
		//FixConstant* getCoeff(int i);


		/** A wrapper for Sollya guessdegree
		 */
		static	void guessDegree(sollya_obj_t fS, sollya_obj_t rangeS, sollya_obj_t targetAccuracy, int& degreeInf, int& degreeSup);

	private:
		/** build an approximation of a certain degree, LSB being already defined, then computes the approx error.
				Essentially a wrapper for Sollya fpminimax() followed by supnorm()
		*/
		void buildApproxFromDegreeAndLSBs();

		/** Build coeff, the vector of coefficients, out of polynomialS, the sollya polynomial
		 	 Constructor code, factored out
		 * */
		void buildFixFormatVector();

        int degree;
    
		sollya_obj_t polynomialS;         /**< The polynomial approximating it */

		sollya_obj_t fixedS;              /**< a constant sollya_obj_t, which indicates that fixed-point formats should be used for fpminimax */
		sollya_obj_t absoluteS;           /**< a constant sollya_obj_t, which indicates that fpminimax must try to minimize the absolute error */

	};

}

#endif
