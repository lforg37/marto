#ifndef KULISCH_DIM_HPP
#define KULISCH_DIM_HPP

#include <array>
#include <cstdint>
#include <type_traits>

#include "tools/static_math.hpp"
#include "helpers/splitting.hpp"
#include "ieeefloats/ieeetype.hpp"

using namespace std;

using hint::Static_Val;
using hint::Static_Ceil_Div;
using hint::ArraySplitter;

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
class FPProd
{
	private:

		using _dim = IEEEDim<WE, WF>;
		constexpr static unsigned int _total_width = _dim::W_Prod;
		constexpr static unsigned int _exp_width = _dim::WE_Prod;
		constexpr static unsigned int _frac_width = _dim::WFF_Prod;

		using  _sign_t = Wrapper<1, false>;
		using _exp_t = Wrapper<_exp_width, false>;
		using _frac_t = Wrapper<_frac_width, false>;
		using _flag_t = Wrapper<1, false>;

		_sign_t _sign;
		_exp_t _exp;
		_frac_t _frac;
		_flag_t _isNaN;
		_flag_t _isInf;

	public:
		FPProd(
				_sign_t const & mult_s,
				_exp_t const & mult_e,
				_frac_t const & mult_m,
				_flag_t const & mult_isNaN,
				_flag_t const & mult_isInf
				) : _sign{mult_s}, _exp{mult_e}, _frac{mult_m}, _isNaN{mult_isNaN}, _isInf{mult_isInf}
		{}

		FPProd(Wrapper<_total_width, false> const & val)
		{
			_sign = val.template get<_total_width - 1>();
			_exp = val.template slice<_total_width - 2, _frac_width>();
			_frac = val.template slice<_frac_width - 1, 0>();
		}

		inline _frac_t getSignificand() const
		{
			return _frac;
		}


		inline _sign_t getSignBit() const
		{
			return _sign;
		}

		inline _exp_t getExp() const
		{
			return _exp;
		}

		inline _flag_t isNaN() const
		{
			return _isNaN;
		}

		inline _flag_t isInf() const
		{
			return _isInf;
		}
};

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
class KulischAcc : public Wrapper<IEEEDim<WE, WF>::ACC_SIZE, false>
{
	private:
		using fpdim = IEEEDim<WE, WF>;
		using storage_type = Wrapper<fpdim::ACC_SIZE, false>;


		inline storage_type const & _storage() const
		{
			return static_cast<storage_type const &>(*this);
		}

	public:
		KulischAcc(storage_type const & in):storage_type{in}
		{}

		Wrapper<1, false> getSignBit() const
		{
			return _storage().template get<fpdim::ACC_SIZE - 1>();
		}

		storage_type const & downcast() const
		{
			return static_cast<storage_type const &>(*this);
		}
};
#endif
