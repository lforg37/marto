#include "add_sub_quire.hpp"
#include "posit_dim.cpp"

//Instantiation of ad_sub_quire for standard sizes
template Quire<8> add_sub_quire<8>(Quire<8>, PositProd<8>, ap_uint<1>);
template Quire<16> add_sub_quire<16>(Quire<16>, PositProd<16>, ap_uint<1>);
template Quire<32> add_sub_quire<32>(Quire<32>, PositProd<32>, ap_uint<1>);
template Quire<64> add_sub_quire<64>(Quire<64>, PositProd<64>, ap_uint<1>);


