#ifndef EXPRESSION_HPP
#define EXPRESSION_HPP

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

#include <tools/static_math.hpp>

using std::add_const;
using std::size_t;

#include <tools/static_math.hpp>

using hint::Static_Val;

#include "fp_number.hpp"

template<typename op1type, typename op2type, template<unsigned int, bool> class Wrapper, template<vec_width, typename, typename> class RoundingOp>
struct BinaryOp
{
	private:
		template<vec_width targetPrecision>
		using op1RoundDim = typename op1type::template roundedDim<targetPrecision>;
		template<vec_width targetPrecision>
		using op2RoundDim = typename op2type::template roundedDim<targetPrecision>;

		template<vec_width targetPrecision>
		using round_op = RoundingOp<targetPrecision, op1RoundDim<targetPrecision>, op2RoundDim<targetPrecision>>;

	public:
		using arg_storage = std::pair<
		typename add_const<op1type>::type,
		typename add_const<op2type>::type
		> const;

		template <vec_width targetPrecision>
		using roundedDim = typename round_op<targetPrecision>::dim;

		template<vec_width maxInternalPrecision>
		static inline FPNumber<roundedDim<maxInternalPrecision>, Wrapper> computeWithTargetPrecision(arg_storage & args)
		{
			auto op1 = args.first.template computeWithTargetPrecision<maxInternalPrecision>();
			auto op2 = args.second.template computeWithTargetPrecision<maxInternalPrecision>();
			return round_op<maxInternalPrecision>::compute(op1, op2);
		}
};

#include "product.hpp"

template<typename op1type, typename op2type, template <unsigned int, bool> class Wrapper>
using ExprProduct = BinaryOp<op1type, op2type, Wrapper, RoundedFPProd>;

template<typename optype, template<unsigned int, bool> class Wrapper, template<vec_width, typename> class RoundingOp>
struct UnaryOp
{
	private:
		template<vec_width targetPrecision>
		using opRoundDim = typename optype::template roundedDim<targetPrecision>;

		template<vec_width targetPrecision>
		using round_op = RoundingOp<targetPrecision, opRoundDim<targetPrecision>>;

	public:
		using arg_storage = typename add_const<optype>::type;
		template <vec_width targetPrecision>
		using roundedDim = typename round_op<targetPrecision>::dim;

		template<vec_width maxInternalPrecision>
		static inline FPNumber<roundedDim<maxInternalPrecision>, Wrapper> computeWithTargetPrecision(arg_storage & arg)
		{
			auto op1 = arg.template computeWithTargetPrecision<maxInternalPrecision>();
			return round_op<maxInternalPrecision>::compute(op1);
		}

};

template<typename optype, template<unsigned int, bool> class Wrapper>
using RoundExpr = UnaryOp<optype, Wrapper, strictRounderOp>;

template<typename optype, template<unsigned int, bool> class Wrapper>
using OppositeExpr = UnaryOp<optype, Wrapper, OppositeOp>;

template<typename optype, template<unsigned int, bool> class Wrapper>
using IdentityExpr = UnaryOp<optype, Wrapper, IdentityOp>;


template<typename dim, template<unsigned int, bool> class Wrapper>
struct LeafFP {
	public:
		using arg_storage = const FPNumber<dim, Wrapper>;
		template <vec_width targetPrecision>
		using roundedDim = dim;

		template<vec_width maxInternalPrecision>
		static inline FPNumber<dim, Wrapper> computeWithTargetPrecision(arg_storage & arg)
		{
			return arg;
		}
};

template<typename T, template<unsigned int, bool> class Wrapper>
struct Expr : public T::arg_storage {
		using arg_storage = typename T::arg_storage;

		template<vec_width targetPrecision>
		using roundedDim = typename T::template roundedDim<targetPrecision>;

		Expr (arg_storage & init):arg_storage{init}
		{}

		template<size_t targetWF>
		inline FPNumber<roundedDim<targetWF>, Wrapper> computeWithTargetPrecision() const
		{
			return T::template computeWithTargetPrecision<targetWF>(static_cast<arg_storage &>(*this));
		}
};

// TODO : find a way to get rid of these ugly "Wrapper"
template<typename dim, template<unsigned int, bool> class Wrapper>
using FPExpr = Expr<LeafFP<dim, Wrapper>, Wrapper>;

template<typename T1, typename T2, template<unsigned int, bool> class Wrapper>
inline Expr<ExprProduct<Expr<T1, Wrapper>, Expr<T2, Wrapper>, Wrapper>, Wrapper> operator*(Expr<T1, Wrapper> const & op1, Expr<T2, Wrapper> const & op2)
{
	return {{op1, op2}};
}


#endif // EXPRESSION_HPP
