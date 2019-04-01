#include "posit_decoder.hpp"

/***
 *
 *  PositEncoding cast
 *
 */
template <int N, int WES>
PositEncoding<N, WES>::operator PositIntermediateFormat<N, WES>() const
{
    #pragma HLS INLINE
    return posit_decoder(*this);
}

template <int N, int WES>
PositEncoding<N, WES>::operator PositProd<N, WES>() const
{
    #pragma HLS INLINE
    return static_cast<PositProd<N, WES> >(static_cast<PositIntermediateFormat<N, WES >>(*this));
}

#include "posit_add.hpp"

template<int N, int WES>
PositEncoding<N, WES> operator+(
        PositEncoding<N, WES> const & lhs,
        PositEncoding<N, WES> const & rhs
    )
{
    #pragma HLS INLINE
    auto lhs_val = static_cast<PositIntermediateFormat<N, WES> >(lhs);
    auto rhs_val = static_cast<PositIntermediateFormat<N, WES> >(rhs);
    return static_cast<PositEncoding<N, WES> >(posit_add_optimized(lhs_val, rhs_val));
}

template<int N, int WES>
PositEncoding<N, WES> operator-(
        PositEncoding<N, WES> const & lhs,
        PositEncoding<N, WES> const & rhs
    )
{
    #pragma HLS INLINE
    auto lhs_val = static_cast<PositIntermediateFormat<N, WES> >(lhs);
    auto rhs_val = static_cast<PositIntermediateFormat<N, WES> >(rhs);
    return static_cast<PositEncoding<N, WES> >(posit_add_optimized(lhs_val, rhs_val, 1));
}

#include "posit_mul.hpp"
template<int N, int WES>
PositProd<N, WES> operator*(
        PositEncoding<N, WES> const & lhs,
        PositEncoding<N, WES> const & rhs
    )
{
    #pragma HLS INLINE
    auto lhs_val = static_cast<PositIntermediateFormat<N, WES> >(lhs);
    auto rhs_val = static_cast<PositIntermediateFormat<N, WES> >(rhs);
    return posit_mul(lhs_val, rhs_val);
}

/***
 *
 * PositIntermediatFormat Conversion
 *
 */

#include "value_prod_conversions.hpp"

template <int N, int WES>
PositIntermediateFormat<N, WES>::operator PositProd<N, WES>() const
{
   #pragma HLS INLINE
   return PositIF_to_PositProd(*this);
}

#include "posit_encoder.hpp"

template<int N, int WES>
PositIntermediateFormat<N, WES>::operator PositEncoding<N, WES>() const
{
    #pragma HLS INLINE
    return posit_encoder(*this);
}

/***
 *
 * PositProd conversion function
 *
 */
template<int N, int WES>
PositProd<N, WES>::operator PositIntermediateFormat<N, WES>() const
{
    #pragma HLS INLINE
    return PositProd_to_PositIF(*this);
}

template<int N, int WES>
PositProd<N, WES>::operator PositEncoding<N, WES>() const
{
    #pragma HLS INLINE
    return PositEncoding<N, WES>{PositProd_to_PositIF(*this)};
}

/***
 *
 * Quire Conversion and arithmetic
 *
 */

#include "quire_to_posit.hpp"
template<int N, int WES, int NB_CARRY>
Quire<N, WES, NB_CARRY>::operator PositIntermediateFormat<N, WES>() const
{
    #pragma HLS INLINE
    return quire_to_posit(*this);
}

template<int N, int WES, int NB_CARRY>
Quire<N, WES, NB_CARRY> operator+(
        Quire<N, WES, NB_CARRY> const & lhs,
        PositProd<N, WES> const & rhs
    )
{
    #pragma HLS INLINE
    return add_sub_quire(lhs, rhs, 0);
}

template<int N, int WES, int NB_CARRY>
Quire<N, WES, NB_CARRY> operator-(
        Quire<N, WES, NB_CARRY> const & lhs,
        PositProd<N, WES> const & rhs
    )
{
    #pragma HLS INLINE
    return add_sub_quire(lhs, rhs, 1);
}

/***
 *
 * SegmentedQuire conversion and arithmetic
 *
 */

template<int N, int WES, int NB_CARRY, int banksize>
SegmentedQuire<N, WES, NB_CARRY, banksize>::operator PositIntermediateFormat<N, WES>() const
{
#pragma HLS INLINE
    auto propagated = propagateCarries(*this);
    return static_cast<PositIntermediateFormat<N, WES> >(propagated);
}


template<int N, int WES, int NB_CARRY, int bankSize>
SegmentedQuire<N, WES, NB_CARRY, bankSize> operator+(
        SegmentedQuire<N, WES, NB_CARRY, bankSize> const & lhs,
        PositProd<N, WES> const & rhs
    )
{
    #pragma HLS INLINE
    return segmented_add_sub_quire(lhs, rhs, 0);
}

template<int N, int WES, int NB_CARRY, int bankSize>
SegmentedQuire<N, WES, NB_CARRY, bankSize> operator-(
        SegmentedQuire<N, WES, NB_CARRY, bankSize> const & lhs,
        PositProd<N, WES> const & rhs
    )
{
    #pragma HLS INLINE
    return segmented_add_sub_quire(lhs, rhs, 1);
}



