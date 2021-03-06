function [maxValue] = TDEmax(MICS,C)

// Output variables initialisation (not found in input variables)
maxValue=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


//Get the maximum possible value of the TDE given the microphones'' position
// 
// USAGE: maxValues = maxTDEs(MICS[,C])
// 
// PARAMETERS:
//   MICS ~ Positions of the microphones (one per row)
//   C ~ Sound speed
// 
//   see also DMIC

// Copyright 2012, Xavier Alameda-Pineda
// INRIA Grenoble Rhône-Alpes
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
// 

//%% Input check
if %nargin<1 then
  error("Usage: maxTDE(MICS[, C])");
end;
if %nargin<2 then
  C = 343.2;
end;

// Variables
NMics = size(MICS,1);
NTDEs = (NMics*(NMics-1))/2;

//%% Computation
maxValue = zeros(NTDEs,1);
for mic1 = 1:NMics
  // !! L.48: Unknown function DMIC not converted, original calling sequence used.
  partialDMIC = DMIC(MICS,MICS(mic1,:));
  for mic2 = mic1+1:NMics
    // Compute the index
    index = TDEis(mic1,mic2,NMics);
    // Compute the maxTDE
    maxValue = mtlb_i(maxValue,index,mtlb_e(partialDMIC,mic2)/C);
  end;
end;

return;
endfunction
