#ifndef IEEE_ADDER_HPP
#define IEEE_ADDER_HPP

#include "ieeefloats/ieeetype.hpp"
#include "ieeefloats/ieee_rounding.hpp"
#include "tools/static_math.hpp"
#include "primitives/lzoc_shifter.hpp"
#include "primitives/lzoc.hpp"
#include "primitives/shifter_sticky.hpp"

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
inline IEEENumber<WE, WF, Wrapper> ieee_add_sub_impl(
	   IEEENumber<WE, WF, Wrapper> in0,
	   IEEENumber<WE, WF, Wrapper> in1,
	   IEEERoundingMode roundingMode = IEEERoundingMode::RoundNearestTieEven
	)
{
	auto exp0 = in0.getExponent();
	auto exp1 = in1.getExponent();

	auto sign0 = in0.getSign();
	auto sign1 = in1.getSign();

	auto frac0 = in0.getFractionnalPart();
	auto frac1 = in1.getFractionnalPart();

	auto expfrac0 = exp0.concatenate(frac0);
	//cerr << "expfrac0 : " << expfrac0.unravel().get_str(2) << endl;
	auto expfrac1 = exp1.concatenate(frac1);
	//cerr << "expfrac1 : " << expfrac1.unravel().get_str(2) << endl;

	auto diff0 = exp0.modularSub(exp1);
	auto diff1 = exp1.modularSub(exp0);

	// Sorting according to exp
	auto swap = expfrac1 > expfrac0;
	//cerr << "swap : " << swap.unravel().get_str(2) << endl;

	auto maxExp = Wrapper<WE, false>::mux(swap, exp1, exp0);
	auto minExp = Wrapper<WE, false>::mux(swap, exp0, exp1);
	auto expdiff = Wrapper<WE, false>::mux(swap, diff1, diff0);

	//cerr << "Expdiff : " << expdiff.unravel().get_str(2) << endl;

	auto effsub = sign0 ^ sign1;

	auto maxSign = Wrapper<1, false>::mux(swap, sign1, sign0);
	auto minSign = Wrapper<1, false>::mux(swap, sign0, sign1);

	auto maxFrac = Wrapper<WF, false>::mux(swap, frac1, frac0);
	auto minFrac = Wrapper<WF, false>::mux(swap, frac0, frac1);

	// Special case detection
	auto maxExpIsZero = maxExp.or_reduction().invert();
	auto minExpIsZero = minExp.or_reduction().invert();
	auto maxExpAllOne = maxExp.and_reduction();
	auto minExpAllOne = minExp.and_reduction();
	auto maxFracIsZero = maxFrac.or_reduction().invert();
	auto minFracIsZero = minFrac.or_reduction().invert();

	//Special case logic
	auto maxIsInfinity = maxExpAllOne & maxFracIsZero;
	auto maxIsNaN = maxExpAllOne & maxFracIsZero.invert();
	auto maxIsZero = maxExpIsZero & maxFracIsZero;
	auto minIsInfinity = minExpAllOne & minFracIsZero;
	auto minIsNaN = minExpAllOne & minFracIsZero.invert();
	auto minIsZero = minExpIsZero & minFracIsZero;

	auto bothSubNormals = maxExpIsZero;
	//cerr << "Both subnormals : " << bothSubNormals.unravel().get_str(2) << endl;
	auto maxIsNormal = maxExpIsZero.invert();
	auto minIsNormal = minExpIsZero.invert();

	auto infinitySub = maxIsInfinity & minIsInfinity & effsub;
	auto resultIsNan = maxIsNaN | minIsNaN |infinitySub;
	auto onlyOneSubnormal = minExpIsZero & maxExpIsZero.invert();
	//cerr << "Only One  subnormal :" << onlyOneSubnormal.unravel().get_str(2) << endl;

	// Reconstruct full fraction
	auto explicitedMaxFrac = maxIsNormal.concatenate(maxFrac);
	auto explicitedMinFrac = minIsNormal.concatenate(minFrac);

	//cerr << "explicit fracs : \nMax : " << explicitedMaxFrac.unravel().get_str(2) << endl <<
	//		"Min : " << explicitedMinFrac.unravel().get_str(2) << endl;

	//alignment
	auto maxShiftVal = Wrapper<WE, false>{WF+3};
	auto allShiftedOut = expdiff > maxShiftVal;
	//cerr << "All shifted out : " << allShiftedOut.unravel().get_str(2) << endl;

	auto shiftValue = Wrapper<WE, false>::mux(
					allShiftedOut,
					maxShiftVal,
					expdiff
				).modularSub(onlyOneSubnormal.template leftpad<WE>());
	//cerr << "Shift value : " << shiftValue.unravel().get_str(2) << endl;

	Wrapper<WF+3, false> extendedMinFrac = explicitedMinFrac.concatenate(Wrapper<2, false>{0});
	//cerr << "extminfrac : " << extendedMinFrac.unravel().get_str(2) << endl;

	auto shiftedMinFracSticky = shifter_sticky(extendedMinFrac, shiftValue);
	auto beforeComp = Wrapper<1, false>{0}.concatenate(shiftedMinFracSticky.template slice<WF + 3, 1>());
	//cerr << "beforeComp : " << beforeComp.unravel().get_str(2) << endl;
	auto shiftedMinFrac = beforeComp ^ Wrapper<WF+4, false>::generateSequence(effsub);
	//cerr << "shiftedMinFrac : " << shiftedMinFrac.unravel().get_str(2) << endl;

	auto stickyMinFrac = shiftedMinFracSticky.template get<0>();
	//cerr << "stickyMinFrac : " << stickyMinFrac.unravel().get_str(2) << endl;

	// Addition
	auto carryIn = effsub & stickyMinFrac.invert();
	//cerr << "carryIn : " << carryIn.unravel().get_str(2) << endl;
	auto extendedMaxFrac = explicitedMaxFrac.concatenate(Wrapper<1, false>{0}).concatenate(carryIn).template leftpad<WF+4>();
	//cerr << "extMaxfrac : " << extendedMaxFrac.unravel().get_str(2) << endl;

	auto signifcandResult = extendedMaxFrac + shiftedMinFrac;
	//cerr << "Signif result : " << signifcandResult.unravel().get_str(2) << endl;
	// Renormalization

	auto isNeg = signifcandResult.template get<WF + 4>();
	auto z1 = signifcandResult.template get<WF+3>();
	auto z0 = signifcandResult.template get<WF+2>();
	//cerr << "Z1 :" << z1.unravel().get_str(2) << endl
	//	<< "Z0  :" << z0.unravel().get_str(2) << endl;

	auto lzcInput = signifcandResult.template slice<WF+3, 1>();
	auto lzc = lzoc_wrapper(lzcInput, {0});

	//cerr << "LZC Input :" << lzcInput.unravel().get_str(2) << endl
	//	<< "LZC  :" << lzc.unravel().get_str(2) << endl;
	constexpr unsigned int lzcsize = hint::Static_Val<WF+3>::_storage;

	static_assert (lzcsize<=WE, "The adder works only for wE > log2(WF).\nAre you sure you need subnormals ?\nIf yes, contact us with your use case, we will be happy to make it work for you.");
	auto subnormalLimitVal = Wrapper<lzcsize, false>{WF+3};

	auto maxExpIsOne = (maxExp == Wrapper<WE, false>{1});

	auto lzcGreaterEqExp = (lzc.template leftpad<WE>() >= maxExp);
	auto lzcSmallerEqExp = (lzc.template leftpad<WE>() <= maxExp);
	auto lzcSmallerMaxVal = lzc < subnormalLimitVal;
	auto fullCancellation = lzcSmallerMaxVal.invert();
	//cerr << "Fulll cancel : " << fullCancellation.unravel().get_str(2) << endl;

	//cerr << "lzcsmaller" << lzcSmallerMaxVal.unravel().get_str(2) << endl;

	auto normalOverflow = z1;
	//cerr << "Normal overflow : " << normalOverflow.unravel().get_str(2) << endl;
	auto lzcOne = z1.invert() & z0;
	auto subnormalOverflow = lzcOne & maxExpIsZero;
	//cerr << "Subnormal overflow : " << subnormalOverflow.unravel().get_str(2) << endl;
	auto cancellation = z1.invert() & z0.invert();

	auto overflow = normalOverflow | subnormalOverflow;

	auto isLeftShiftLZC = overflow |
			(lzcSmallerMaxVal.invert() & bothSubNormals) |
			(cancellation & maxIsNormal & lzcSmallerEqExp) |
			(maxIsNormal & lzcSmallerMaxVal.invert());
	auto isLeftShiftExp = lzcSmallerMaxVal & lzcGreaterEqExp & maxIsNormal;

	auto shiftFirstStage = Wrapper<lzcsize, false>::mux(isLeftShiftLZC, lzc, Wrapper<lzcsize, false>{1});
	auto normalisationShiftVal = Wrapper<lzcsize, false>::mux(isLeftShiftExp, maxExp.template slice<lzcsize-1, 0>(), shiftFirstStage);

	//cerr << "Normalisation shift : " << normalisationShiftVal.unravel().get_str(2) << endl;
	auto normalisedSignif = signifcandResult << normalisationShiftVal;
	//cerr << "Normalised signif : " << normalisedSignif.unravel().get_str(2) << endl;
	auto significandPreRound = normalisedSignif.template slice<WF+2, 3>();
	//cerr << "Preround signif : " << significandPreRound.unravel().get_str(2) << endl;
	auto lsb = normalisedSignif.template get<3>();
	//cerr << "lsb : " << lsb.unravel().get_str(2) << endl;
	auto roundBit = normalisedSignif.template get<2>();
	//cerr << "roundBit : " << roundBit.unravel().get_str(2) << endl;
	auto sticky = stickyMinFrac | normalisedSignif.template slice <1, 0>().or_reduction();
	//cerr << "sticky : " << sticky.unravel().get_str(2) << endl;

	auto deltaExpIsZero = z1.invert() & (z0 ^ bothSubNormals);
	//cerr << "DeltaExpIsZero : " << deltaExpIsZero.unravel().get_str(2) << endl;
	auto deltaExpIsMinusOne = z1 | (z0 & bothSubNormals);
	//cerr << "deltaExpIsMinusOne : " << deltaExpIsMinusOne.unravel().get_str(2) << endl;
	auto deltaExpIsLZC = (z1 | z0 | bothSubNormals).invert() & lzcSmallerEqExp & lzcSmallerMaxVal;
	//cerr << "deltaExpIsLZC : " << deltaExpIsLZC.unravel().get_str(2) << endl;
	auto deltaExpExp = (deltaExpIsLZC | deltaExpIsZero | deltaExpIsMinusOne).invert();

	auto deltaExpCin = deltaExpExp | deltaExpIsMinusOne | deltaExpIsLZC;
	auto deltaBigPartIsZero = deltaExpIsZero | deltaExpIsMinusOne;

	auto deltaExpUnmasked = Wrapper<WE, false>::mux(
				deltaExpIsLZC,
				(lzc.modularSub(Wrapper<lzcsize, false>{1})).template leftpad<WE>().invert(),
				maxExp.invert()
			);

	auto maskSequence = Wrapper<WE, false>::generateSequence(deltaBigPartIsZero.invert());
	auto deltaExpBeforeCorrection = deltaExpUnmasked & maskSequence;

	//cerr << "DeltaExp : " << deltaExpBeforeCorrection.unravel().get_str(2) << endl;
	//cerr << "DeltaEXP cin : " << deltaExpCin.unravel().get_str(2) <<endl;

	auto expPreRound = maxExp.addWithCarry(deltaExpBeforeCorrection, deltaExpCin).template slice<WE-1, 0>();
	//cerr << "expPreRound : " << expPreRound.unravel().get_str(2) << endl;
	auto expSigPreRound = expPreRound.concatenate(significandPreRound);
	//cerr << "expSigPreRound : " << expSigPreRound.unravel().get_str(2) << endl;

	auto roundUpBit = ieee_getRoundBit(maxSign, lsb, roundBit, sticky, roundingMode);
	//auto roundUpBit = roundBit & (sticky | lsb);
	auto expSigRounded = expSigPreRound.modularAdd(roundUpBit.template leftpad<WE+WF>());
	auto finalExp = expSigRounded.template slice<WF+WE-1, WF>();

	auto resultIsZero = fullCancellation & finalExp.or_reduction().invert();
	//cerr << "ResIsZero : " << resultIsZero.unravel().get_str(2) << endl;
	auto resultIsInf = resultIsNan.invert() & (
					(maxIsInfinity & minIsInfinity & effsub.invert()) |
					(maxIsInfinity ^ minIsInfinity) |
					finalExp.and_reduction()
				);

	auto constInfNan = Wrapper<WE, false>::generateSequence({1}).concatenate(Wrapper<WF, false>::generateSequence(resultIsNan));
	auto finalRes = Wrapper<WE+WF, false>::mux(resultIsNan | resultIsInf, constInfNan, expSigRounded);

	auto bothZeros = maxIsZero & minIsZero;
	auto signBothZero = minSign & maxSign;

	auto signR = resultIsNan.invert() & // NaN forces sign to be zero
			(resultIsZero & signBothZero.invert()).invert() & // If summing two zeros of opposite sign, set result to zero
			maxSign;

	return {signR.concatenate(finalRes)};
}
#endif // IEEE_ADDER_HPP
