#ifndef SHIFTER_STICKY_HPP
#define SHIFTER_STICKY_HPP

#include <iostream>
#include <type_traits>

#include "ap_int.h"

#include "shifter.hpp"

template<int N, int S>
//N : Power of 2 of the size of the whole LZOC,
//S : size of the shift to consider
inline ap_uint<(1 << N) + 1> sticky_shifter_stage(
                ap_uint<(1<<N) + 1> input,
                ap_uint<S> count,
                ap_uint<1> fill_bit = 0,
                typename std::enable_if<ShifterStageInfo<S>::NeedsRecursion>::type* = 0
        )
{
        #pragma HLS INLINE
        ap_int<1<<(S-1)> padding_s = static_cast<ap_int<1> >(fill_bit);
        ap_uint<1<<(S-1)> padding = padding_s;
        ap_uint<1> stageNeedsShift = count[S-1];
        ap_uint<S-1> countnext = count.range(S-2, 0);

        ap_uint<1> sticky_in = input[0];
        ap_uint<(1 << N) - (1 << (S-1))> low = input.range((1 << N) - (1 << (S-1)), 1);
        ap_uint<(1 << (S-1))> high = input.range((1 << N), 1<< (S-1) + 1);
        ap_uint<1<<N> next_stage_input;
        ap_uint<1> cur_sticky;

        if (stageNeedsShift) {
                cur_sticky = sticky_in or high.or_reduce();
                next_stage_input = low.concat(padding);
        } else {
                cur_sticky = sticky_in;
                next_stage_input = input.range(1 << N, 1);
        }
        return shifter_stage<N, S-1>(next_stage_input.concat(cur_sticky), countnext, fill_bit);
}

template<int N, int S>
inline ap_uint<(1 << N) + 1> shifter_sticky_stage(
                ap_uint<(1<<N) + 1> input,
                ap_uint<S> count,
                ap_uint<1> fill_bit = 0,
                typename std::enable_if<ShifterStageInfo<S>::IsFinalStage>::type* dummy = 0
        )
{
        #pragma HLS INLINE
        ap_uint<(1<<N) + 1> result;
        if (count[0] == 1) {
                ap_uint<(1<<N) - 1> low = input.range((1<<N) - 1, 1);
                ap_uint<1> sticky_out = input[0] or input[1 << N];
                result = static_cast<ap_uint<1<<N> >(low.concat(fill_bit)).concat(sticky_out);
        } else {
                result = input;
        }
        return result;
}

template<int N, bool isRightShift>
ap_uint<(1<<N)+1> shifter_sticky(
                ap_uint<1<<N> input,
                ap_uint<N-1> count,
                ap_uint<1> fill_bit = 0)
{
        #pragma HLS INLINE
        ap_uint<1<<N> fin_input;
        if (isRightShift) {
            fin_input = input.reverse();
        } else {
            fin_input = input;
        }
        ap_uint<(1<<N) + 1> shiftstick = shifter_sticky_stage<N, N-1>(fin_input.concat(ap_uint<1>{0}), count, fill_bit);
        ap_uint<(1<<N) + 1> ret;
        if (isRightShift) {
            ap_uint<1<<N> high = shiftstick.range(1<<N, 1);
            ap_uint<1<<N> reverse = high.reverse();
            ap_uint<1> sticky = shiftstick[0];
            ret = reverse.concat(sticky);
        } else {
            ret = shiftstick;
        }
        return ret;
}

#endif // SHIFTER_STICKY_HPP
