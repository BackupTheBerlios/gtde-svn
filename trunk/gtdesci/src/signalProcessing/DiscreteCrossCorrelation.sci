function [%corr] = DiscreteCrossCorrelation(signal1,signal2,maxLag)

// Output variables initialisation (not found in input variables)
%corr=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


//Discrete cross-correlation
// 
// USAGE: Cross-Correlation = DiscreteCrossCorrelation(FirstSignal,SecondSignal,MaxLag)
// 
// PARAMETERS:
//   FirstSignal ~ values of the first signal
//   SecondSignal ~ values of the second signal
//   MaxLag ~ the maximum value of the lag we want to compute
// 
// RETURN VALUE:
//   Cross-Correlation ~ values of the cross-correlation signal
// 
// DESCRIPTION:
//   This function computes the discrete cross-corrlation of two 1D signals.
// 

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

//%% Check input
if %nargin<3 then
  error("Usage: Cross-Correlation = DiscreteCrossCorrelation(FirstSignal,SecondSignal,MaxLag)");
end;

//%% Initialize output
corr = zeros(1,mtlb_a(2*maxLag,1));
// How many samples we use to compute each of the cross-correlation
// values
chunkSize = mtlb_s(max(size(signal1)),maxLag);

//%% Computation
// Positive values of the delay
for xcPos = mtlb_imp(2,mtlb_a(maxLag,1))
  %corr = mtlb_i(%corr,mtlb_a(maxLag,xcPos),mtlb_sum(mtlb_e(signal1,mtlb_imp(xcPos,mtlb_a(mtlb_s(chunkSize,1),xcPos))) .*mtlb_e(signal2,mtlb_imp(1,chunkSize))));
end;
// Negative values of the delay
for xcPos = mtlb_imp(2,mtlb_a(maxLag,1))
  %corr = mtlb_i(%corr,mtlb_s(mtlb_a(maxLag,2),xcPos),mtlb_sum(mtlb_e(signal2,mtlb_imp(xcPos,mtlb_a(mtlb_s(chunkSize,1),xcPos))) .*mtlb_e(signal1,mtlb_imp(1,chunkSize))));
end;
// At 0
%corr = mtlb_i(%corr,mtlb_a(maxLag,1),mtlb_sum(mtlb_e(signal1,mtlb_imp(1,chunkSize)) .*mtlb_e(signal2,mtlb_imp(1,chunkSize))));

return;

endfunction
