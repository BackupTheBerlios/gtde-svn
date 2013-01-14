function [CrossCorrelation,Derivative,Curvature] = InterpolateXCorrelation(FirstSignal,SecondSignal,Delays,SamplingPeriod,MaxDelay)

// Output variables initialisation (not found in input variables)
CrossCorrelation=[];
Derivative=[];
Curvature=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


//InterpolateXCorrelation interpolates the cross-correlation function 
//of two discrete signals
// 
// USAGE: [CrossCorrelation Derivative Curvature] = 
//   InterpolateXCorrelation(FirstSignal,SecondSignal,Delays,SamplingPeriod)
// 
// PARAMETERS:
//   FirstSignal ~ values of the first signal
//   SecondSignal ~ values of the second signal
//   Delays ~ delays in which we want to evaluate the cross-correlation function
//   SamplingPeriod ~ sampling period of the signals
// 
// RETURN VALUE:
//   CrossCorrelation ~ values of the cross-correlation function at delays
//   Derivatives ~ values of its derivative at delays
//   Curvature ~ values of its second derivative at delays
// 
// DESCRIPTION:
//   This function interpolates the cross-correlation function at times Delays
//   of the signals FirstSignal and SecondSignal (sampled at SamplingPeriod). 
//   It does that assuming some polynomial interpolation at each time interval.
// 
//   This function will then compute:
//   1) The sequence of coefficients of each polynomial interpolation.
//   2) The cross-correlation of these sequences of coefficients.
//   3) The interpolation of the cross-correlation function at each of the Delays.
// 
// REFERENCES:
//     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
//     estimation-based sound source localisation (gTDESSL). Research Report 
//     RR-7988, INRIA, June 2012.
// 
//   see also PolynomialInterpolationCoefficients, PolynomialCoefficientsCrossCorrelation and performCCInterpolation

// Copyright 2012, Xavier Alameda-Pineda
// INRIA Grenoble Rhone-Alpes
// E-mail: xavi.alameda@gmail.com
// 
// This is part of the gtde program.
// 
// gtde is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

getd('../src')

//%% Input check
if %nargin<4 then
  error("Usage: [CrossCorrelation Derivative] = InterpolateXCorrelation(FirstSignal,SecondSignal,Delays,SamplingPeriod[,MaxDelay]");
end;

//%% Compute the interpolation polynomial coefficients of the signals

FirstPolynomialCoefficie = PolyInterpolation(FirstSignal,SamplingPeriod);

SecondPolynomialCoeffici = PolyInterpolation(SecondSignal,SamplingPeriod);

//%% Compute the cross-correlation of the coefficients
if %nargin<5 then

  PCCC = PolyCoefXCorrelation(FirstPolynomialCoefficie,SecondPolynomialCoeffici);
else
  MaxLag = ceil(MaxDelay/SamplingPeriod);


  PCCC = PolyCoefXCorrelation(FirstPolynomialCoefficie,SecondPolynomialCoeffici,MaxLag);
end;


//%% Interpolate the cross-correlation at """"Delays""""
%v0 = size(Delays);CrossCorrelation = zeros(%v0(1),%v0(2));
%v0 = size(Delays);Derivative = zeros(%v0(1),%v0(2));
%v0 = size(Delays);Curvature = zeros(%v0(1),%v0(2));

// !! L.84: Matlab function numel not yet converted, original calling sequence used.

for dd = mtlb_imp(1,length(Delays))


  [%v0,%v1,%v2] = XCorrInterpolation(PCCC,mtlb_e(Delays,dd),SamplingPeriod);
  CrossCorrelation = mtlb_i(CrossCorrelation,dd,%v0);
  Derivative = mtlb_i(Derivative,dd,%v1);
  Curvature = mtlb_i(Curvature,dd,%v2);
end;

return;

endfunction
