#include <ap_int.h>
#include "hint.hpp"
#include "ieeefloats/ieeetype.hpp"
#include "ieeefloats/ieee_adder.hpp"

constexpr unsigned int WE = IEEE_WE;
constexpr unsigned int WF = IEEE_WF;
constexpr unsigned int ENC_SIZE = WE + WF + 1;

ap_uint<ENC_SIZE> ieee_comp_add(
	IEEENumber<WE, WF, hint::VivadoWrapper> in1,
	IEEENumber<WE, WF, hint::VivadoWrapper> in2
)
{
	auto add = ieee_add_sub_impl(in1, in2);
	return add.unravel();
}

int main(void)
{
	ieee_comp_add({{17}}, {{78}});
	return 0;
}
