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

constexpr int ceil2Power(int N)
{
	if(N==0){
		return 0;
	}
	else{
		int tmp_N = N;
		int result = 1;
		while(tmp_N!=0){
			tmp_N = tmp_N>>1;
			result = result<<1;
		}
		if((result>>1) == N){
			return N;
		}
		else{
			return result;	
		}
	}
}

constexpr int ceilLog2(int N)
{
	return (int) ceil(log2(N));
}

template<int N>
class PositDim {
	public:
	static constexpr int WES = get2Power(N>>3); 
	static constexpr int WE = (N-2)*(1 << WES); 
	static constexpr int WF = N - (WES+3);
	//Quire dimension
	static constexpr int WQ = (N*N) >> 1;

	static constexpr int ValSize = 3 + WE + WF;
	static constexpr int ExtQuireSize = WQ + 1;
	static constexpr int ProdSignificandSize = 2*WF + 2;
	static constexpr int ProdExpSize = WE + 1;
	static constexpr int ProdSize = 1 + ProdExpSize + ProdSignificandSize;

	static constexpr bool HAS_ES = (WES > 0);
};

//One bit is NaR + quire
template<int N>
class Quire
{
	//Storage :
	// isNar Sign Carry 2sCompValue
	public:
		Quire(ap_uint<PositDim<N>::ExtQuireSize> val):_val(val){}

		ap_uint<PositDim<N>::WQ-N> getQuireValue()
		{
			return _val.range(PositDim<N>::WQ-N -1, 0);
		}

		ap_uint<N-1> getCarry()
		{
			return _val.range(PositDim<N>::WQ-1 -1, 
					PositDim<N>::WQ-N);
		}

		ap_uint<PositDim<N>::ExtQuireSize-2> getQuireWithoutNaR()
		{
			return _val.range(PositDim<N>::ExtQuireSize-1 -1, 0);
		}

		ap_uint<1> getSignBit()
		{
			return _val[PositDim<N>::ExtQuireSize-1 -1];
		}


		ap_uint<1> getIsNaR()
		{
			return _val[PositDim<N>::ExtQuireSize -1];
		}

		ap_uint<1> getBit(unsigned int i)
		{
			return _val[i];
		}

	private:
		ap_uint<PositDim<N>::ExtQuireSize> _val;	
};


// One bit isNar + WE+1 + 2(WF+1) 
template<int N>
class PositProd
{
	//Storage :
	// isNar Exp Sign ImplicitBit Fraction
	public:
		PositProd(ap_uint<PositDim<N>::ProdSize> val):_val(val){}

		ap_uint<2*PositDim<N>::WF+2> getSignificand()
		{
			return _val.range(2*PositDim<N>::WF+2 -1, 0);
		}

		ap_uint<1> getSignBit()
		{
			return _val[2*PositDim<N>::WF+2];
		}

		ap_uint<PositDim<N>::WE> getExp()
		{
			return _val.range(PositDim<N>::WE+1+2*PositDim<N>::WF+2 -1, 
					2*PositDim<N>::WF+2+1);
		}

		ap_uint<1> getIsNaR()
		{
			return _val[PositDim<N>::WE+1+2*PositDim<N>::WF+2 -1];
		}

		ap_uint<1> getBit(unsigned int i)
		{
			return _val[i];
		}

	private:
		ap_uint<PositDim<N>::ProdSize> _val;	
};

template<int N>
using PositEncoding = ap_uint<N>;

template<int N>
class PositValue
{
	//Storage :
	// isNar Exp Sign ImplicitBit Fraction
	public:
		PositValue(
				ap_uint<1> isNar,
				ap_uint<PositDim<N>::WE> exp,
				ap_uint<1> sign,
				ap_uint<1> implicit_bit,
				ap_uint<PositDim<N>::WF> fraction)
		{
			ap_uint<1+PositDim<N>::WE> tmp = isNar.concat(exp); 
			ap_uint<2> frac_lead = sign.concat(implicit_bit);
			ap_uint<2+PositDim<N>::WF> full_frac = frac_lead.concat(fraction);
			_val = tmp.concat(full_frac);
		}
				
		PositValue(ap_uint<PositDim<N>::ValSize> val):_val(val){}

		ap_uint<PositDim<N>::WF+1> getSignificand()
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

		ap_uint<1> getBit(unsigned int i)
		{
			return _val[i];
		}

		void printContent(){
			fprintf(stderr, "isNaR: %d\n", (int) this->getIsNaR());
			
			fprintf(stderr, "exp: ");
			for(int i=PositDim<N>::WE -1; i>=0; i--){
				fprintf(stderr, "%d", (int) (this->getExp())[i]);
			}
			fprintf(stderr, "\n");

			fprintf(stderr, "sign: %d\n", (int) this->getSignBit());

			fprintf(stderr, "significand: %d.", (int) (this->getSignificand())[PositDim<N>::WF+1 -1]);
			for(int i=PositDim<N>::WF+1-1 -1; i>=0; i--){
				fprintf(stderr, "%d", (int) (this->getSignificand())[i]);
			}
			fprintf(stderr, "\n");

		}

	private:
		ap_uint<PositDim<N>::ValSize> _val;	
};

#endif
