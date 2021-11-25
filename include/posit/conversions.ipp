#pragma once

#include "posit_decoder.hpp"

/***
 *
 *  PositEncoding cast
 *
 */
template <unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
inline PositEncoding<N, WES, Wrapper>::operator PositIntermediateFormat<N, WES, Wrapper, true>() const
{
	return posit_decoder(*this);
}

template <unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
inline PositEncoding<N, WES, Wrapper>::operator PositProd<N, WES, Wrapper>() const
{
	return static_cast<PositProd<N, WES, Wrapper> >(static_cast<PositIntermediateFormat<N, WES, Wrapper, true>>(*this));
}

#include "posit_add.hpp"

template<unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
inline PositEncoding<N, WES, Wrapper> operator+(
		PositEncoding<N, WES, Wrapper> const lhs,
		PositEncoding<N, WES, Wrapper> const rhs
	)
{
	auto lhs_val = static_cast<PositIntermediateFormat<N, WES, Wrapper, true> >(lhs);
	auto rhs_val = static_cast<PositIntermediateFormat<N, WES, Wrapper, true> >(rhs);
	return static_cast<PositEncoding<N, WES, Wrapper> >(posit_add(lhs_val, rhs_val));
}

template<unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
inline PositEncoding<N, WES, Wrapper> operator-(
		PositEncoding<N, WES, Wrapper> const lhs,
		PositEncoding<N, WES, Wrapper> const rhs
	)
{
	auto lhs_val = static_cast<PositIntermediateFormat<N, WES, Wrapper, true> >(lhs);
	auto rhs_val = static_cast<PositIntermediateFormat<N, WES, Wrapper, true> >(rhs);
	return static_cast<PositEncoding<N, WES, Wrapper> >(posit_add(lhs_val, rhs_val, 1));
}

#include "posit_mul.hpp"
template<unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
inline PositProd<N, WES, Wrapper> operator*(
		PositEncoding<N, WES, Wrapper> const lhs,
		PositEncoding<N, WES, Wrapper> const rhs
	)
{
	auto lhs_val = static_cast<PositIntermediateFormat<N, WES, Wrapper, true> >(lhs);
	auto rhs_val = static_cast<PositIntermediateFormat<N, WES, Wrapper, true> >(rhs);
	return posit_mul(lhs_val, rhs_val);
}

/***
 *
 * PositIntermediatFormat Conversion
 *
 */

#include "value_prod_conversions.hpp"

template <unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, true>::operator PositProd<N, WES, Wrapper>() const
{
   return PositIF_to_PositProd(*this);
}

#include "posit_encoder.hpp"

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, false>::operator PositEncoding<N, WES, Wrapper>() const
{
	return posit_encoder(*this);
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, true>::operator PositEncoding<N, WES, Wrapper>() const
{
	return posit_encoder(*this);
}

/***
 *
 * PositProd conversion function
 *
 */
template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositProd<N, WES, Wrapper>::operator PositIntermediateFormat<N, WES, Wrapper, false>() const
{
	return PositProd_to_PositIF(*this);
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositProd<N, WES, Wrapper>::operator PositEncoding<N, WES, Wrapper>() const
{
	return PositEncoding<N, WES, Wrapper>{PositProd_to_PositIF(*this)};
}


/***
 *
 * Quire Conversion and arithmetic
 *
 */

#include "quire_to_posit.hpp"
template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY>
inline Quire<N, WES, Wrapper, NB_CARRY>::operator PositIntermediateFormat<N, WES, Wrapper, false>() const
{
	return quire_to_posit(*this);
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY>
inline Quire<N, WES, Wrapper, NB_CARRY> operator+(
		Quire<N, WES, Wrapper, NB_CARRY> const lhs,
		PositProd<N, WES, Wrapper> const rhs
	)
{
	return add_sub_quire(lhs, rhs, 0);
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY>
inline Quire<N, WES, Wrapper, NB_CARRY> operator-(
		Quire<N, WES, Wrapper, NB_CARRY> const lhs,
		PositProd<N, WES, Wrapper> const rhs
	)
{
	return add_sub_quire(lhs, rhs, 1);
}

/***
 *
 * SegmentedQuire conversion and arithmetic
 *
 */

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int banksize>
inline SegmentedQuire<N, WES, Wrapper, NB_CARRY, banksize>::operator PositIntermediateFormat<N, WES, Wrapper, false>() const
{
	auto propagated = propagateCarries(*this);
	return static_cast<PositIntermediateFormat<N, WES, Wrapper, false> >(propagated);
}


template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int banksize>
inline SegmentedQuire<N, WES, Wrapper, NB_CARRY, banksize> operator+(
		SegmentedQuire<N, WES, Wrapper, NB_CARRY, banksize> const lhs,
		PositProd<N, WES, Wrapper> const rhs
	)
{
	return segmented_add_sub_quire(lhs, rhs, 0);
}
template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int banksize>
inline SegmentedQuire<N, WES, Wrapper, NB_CARRY, banksize> operator-(
		SegmentedQuire<N, WES, Wrapper, NB_CARRY, banksize> const lhs,
		PositProd<N, WES, Wrapper> const rhs
	)
{
	return segmented_add_sub_quire(lhs, rhs, 1);
}



