#ifndef POSIT_DIM_TPP
#define POSIT_DIM_TPP

using namespace std;

#include "ap_int.h"

//N.B.: We are using int instead of size_t because of ap_uint is templatized
//=====  with int

constexpr int get2Power(int N)
{
	return (N<=1) ? 0 : 1 + get2Power(N >> 1); 
}

template<int N>
class PositDim {
	public:
	static constexpr int WE = get2Power(N>>3); 
	static constexpr int WF = N - (WE+3);
	//Quire dimension
	static constexpr int WQ = (N*N) >> 1;

	static constexpr int ValSize = 3 + WE + WF;
	static constexpr int ExtQuireSize = WQ + 1;
	static constexpr int ProdSize = 4 + WE + 2*WF;

	static constexpr bool HAS_ES = (WE > 0);
};

//One bit is NaR + quire
template<int N>
using Quire = ap_uint<PositDim<N>::ExtQuireSize>;

// One bit isNar + WE+1 + 2(WF+1) 
template<int N>
using PositProd = ap_uint<PositDim<N>::ProdSize>;

template<int N>
using PositEncoding = ap_uint<N>;

template<int N>
class PositValue
{
	//Storage :
	// isNar Exp Sign ImplicitBit Fraction
	public:
		PositValue(ap_uint<PositDim<N>::ValSize> val):_val(val){}

		ap_uint<PositDim<N>::WF+2> getSignificand()
		{
			return _val.range(PositDim<N>::WF, 0);
		}

		ap_uint<1> getSignBit()
		{
			return _val[PositDim<N>::WF + 1];
		}

		ap_uint<PositDim<N>::WE> getExp()
		{
			return _val.range(PositDim<N>::WF + 1 + PositDim<N>::WE, 
					PositDim<N>::WF+2);
		}

		ap_uint<1> getIsNaR()
		{
			return _val[PositDim<N>::WF+2+PositDim<N>::WE];
		}

	private:
		ap_uint<PositDim<N>::ValSize> _val;	
};

#endif
