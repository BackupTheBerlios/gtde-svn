function [cc,ccd,ccdd] = XCorrInterpolation(PCCC,Delay,SamplingPeriod,ZeroIndex)

// Output variables initialisation (not found in input variables)
cc=[];
ccd=[];
ccdd=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


//ccInterpolation Cross-correlation function interpolation
// 
// USAGE: [cc ccd ccdd] = XCorrInterpolation(PCCC,Delay,T[,ZI])
// 
// PARAMETERS:
//  PCCC ~ cross-correlation signals of the polynomial coefficients.
//  Delay ~ delay is the point(s) in which we want to estimate the signals''
//     cross-correlation function.
//  T ~ the signals'' sampling period.
//  ZI ~ the index on the PCCC representing no delay (the default value is ceil(L/2), where L is the
//     size of the coefficients cross-correlation function.)
// 
// RETURN VALUE:
//  cc,ccd,ccdd ~ The value of the cross-correlation function, its first and
//     its second derivative.
// 
// DESCRIPTION:
//     This function estimates the value of the cross-correlation function
//     at Delay from from the polynomial coefficients'' cross-correlation 
//     signals PCCC.
//    
// REFERENCES:
//     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
//     estimation-based sound source localisation (gTDESSL). Research Report 
//     RR-7988, INRIA, June 2012.
// 
//     see also K1, K2, InterpolateCrossCorrelation

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

// Input check
if %nargin<3 then
  error("Usage: [cc ccd ccdd] = XCorrInterpolation(PCCC,Delay,SamplingPeriod[,ZeroIndex])");
end;
if %nargin<4 then
  // L.55: M2SCI found: 3 > size(size(PCCC),2),
  // So result is set to 1.
  ZeroIndex = ceil(1/2);
end;

// Output value
%v0 = size(Delay);cc = zeros(%v0(1),%v0(2));
%v0 = size(Delay);ccd = zeros(%v0(1),%v0(2));
%v0 = size(Delay);ccdd = zeros(%v0(1),%v0(2));

// Index in the cross-correlation
ccIndex = floor(Delay/SamplingPeriod);
// Remaining of the division (0 <= tau < SamplingPeriod)
tau = mtlb_s(Delay,ccIndex*SamplingPeriod);

// Interpolation polynomial degrees
D1 = size(PCCC,1)-1;
D2 = size(PCCC,2)-1;

for d1 = 0:D1
  for d2 = 0:D2
    try 
      // Value of the first interpolating function
      // !! L.76: Unknown function K1 not converted, original calling sequence used.
      [k1,k1d,k1dd] = K1(d1,d2,tau,SamplingPeriod);
      // Value of the second interpolation function
      // !! L.78: Unknown function K2 not converted, original calling sequence used.
      [k2,k2d,k2dd] = K2(d1,d2,tau,SamplingPeriod);
      // First cross-correlation value
      // !! L.80: Matlab function squeeze not yet converted, original calling sequence used.
      R1 = mtlb_t(squeeze(PCCC(d1+1,d2+1,mtlb_a(mtlb_a(ccIndex,1),ZeroIndex))));
      // Second cross-correlation value
      // !! L.82: Matlab function squeeze not yet converted, original calling sequence used.
      R2 = mtlb_t(squeeze(PCCC(d1+1,d2+1,mtlb_a(ccIndex,ZeroIndex))));
      // Cross-correlation''s value
      cc = mtlb_a(mtlb_a(cc,R1 .*k1),R2 .*k2);
      // Cross-correlation''s derivative
      ccd = mtlb_a(mtlb_a(ccd,R1 .*k1d),R2 .*k2d);
      // Cross-correlation''s second derivative
      ccdd = mtlb_a(mtlb_a(ccdd,R1 .*k1dd),R2 .*k2dd);
    catch  exception
      disp(''Exception'');
    end;
  end;
end;
return;

endfunction
