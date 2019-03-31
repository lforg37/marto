#pragma once

#include <iostream>
#include <type_traits>

#include "ap_int.h"

#include "shifter.hpp"
#include "posit_dim.hpp"
#include "utils.hpp"

template<int IS, int S>
//IS : Input Size (including sticky bit),
//S : size of shift counter
ap_uint<IS> shifter_sticky_stage(
                ap_uint<IS> input,
                ap_uint<S> count,
                ap_uint<1> fill_bit = 0,
                typename std::enable_if<ShifterStageInfo<S>::NeedsRecursion>::type* = 0,
                typename std::enable_if<((IS-1) >= (1 << (S-1)))>::type * = 0
        )
{
        #pragma HLS INLINE
        ap_int<1<<(S-1)> padding_s = static_cast<ap_int<1> >(fill_bit);

        ap_uint<1<<(S-1)> padding = padding_s;
        ap_uint<1> stageNeedsShift = count[S-1];
        ap_uint<S-1> countnext = count.range(S-2, 0);

        ap_uint<1> sticky_in = input[0];
        ap_uint<IS - (1 << (S-1))> low = input.range(IS - 1 - (1 << (S-1)), 1);
        ap_uint<(1 << (S-1))> high = input.range(IS - 1 , IS - (1 << (S-1)));
        ap_uint<IS> next_stage_input;
        if (stageNeedsShift) {
                ap_uint<1> cur_sticky = sticky_in or high.or_reduce();
                next_stage_input = static_cast<ap_uint<IS-1> >(low.concat(padding)).concat(cur_sticky);
        } else {
                next_stage_input = input;
        }
        return shifter_sticky_stage<IS, S-1>(next_stage_input, countnext, fill_bit);
}

template<int IS, int S>
ap_uint<IS> shifter_sticky_stage(
                ap_uint<IS> input,
                ap_uint<S> count,
                ap_uint<1> fill_bit = 0,
                typename std::enable_if<ShifterStageInfo<S>::IsFinalStage>::type* = 0,
                typename std::enable_if<((IS-1) >= (1 << (S-1)))>::type * = 0
        )
{
        #pragma HLS INLINE
        ap_uint<IS> result;
        if (count[0] == 1) {
                ap_uint<IS - 2> low = input.range(IS - 2, 1);
                ap_uint<1> sticky_out = input[0] or input[IS-1];
                result = static_cast<ap_uint<IS-1> >(low.concat(fill_bit)).concat(sticky_out);
        } else {
                result = input;
        }
        return result;
}

template<int IS, int S>
inline ap_uint<IS> shifter_sticky_stage(
        ap_uint<IS> input,
        ap_uint<S> count,
        ap_uint<1> fill_bit = 0,

        typename std::enable_if<((IS-1) < (1 << (S-1)))>::type * = 0
    )
{
    #pragma HLS INLINE
    constexpr int nb_null_shift = S - Static_Val<IS-1>::_log2;
    ap_uint<nb_null_shift> shift_weights_will_zero = count.range(S - 1,
                                                                 S - 1 - nb_null_shift);
    ap_uint<S-nb_null_shift> next_count = count.range(S-2-nb_null_shift, 0);

    ap_uint<1> stageNeedsShift = shift_weights_will_zero.or_reduce();
    ap_uint<IS> ret;
    if (stageNeedsShift) {
        ap_uint<1> sticky = input.or_reduce();
        ap_int<IS-1> high = static_cast<ap_int<1> >(fill_bit);
        ret = high.concat(sticky);
    } else {
        ret = shifter_sticky_stage<IS, S-nb_null_shift>(input, next_count, fill_bit);
    }
    return ret;
}

template<int IS, int S, bool isRightShift>
ap_uint<IS+1> shifter_sticky(
                ap_uint<IS> input,
                ap_uint<S> count,
                ap_uint<1> fill_bit = 0
    )
{
        #pragma HLS INLINE
        ap_uint<IS> fin_input;
        if (isRightShift) {
            fin_input = input.reverse();
        } else {
            fin_input = input;
        }
        ap_uint<IS + 1> init_sticky = fin_input.concat(ap_uint<1>{0});
        ap_uint<IS + 1> shiftstick = shifter_sticky_stage<IS+1, S>(init_sticky, count, fill_bit);
        ap_uint<IS + 1> ret;
        if (isRightShift) {
            ap_uint<IS> high = shiftstick.range(IS, 1);
            ap_uint<1> sticky = shiftstick[0];
            ap_uint<IS> reverse = high.reverse();
            ret = reverse.concat(sticky);
        } else {
            ret = shiftstick;
        }
        return ret;
}
