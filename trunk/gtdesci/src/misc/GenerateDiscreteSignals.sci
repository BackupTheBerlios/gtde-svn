function [signals,times] = GenerateDiscreteSignals(Position,MICS,samplingPeriod,Length,signal,C)

// Output variables initialisation (not found in input variables)
signals=[];
times=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


//Generates discrete signals from their analytic expressions
// 
// USAGES: 
//   [signals times] = GenerateDiscreteSignals(Position,MICS,samplingPeriod,Length,signalChoice[,C])
//   [signals times] = GenerateDiscreteSignals(samplingPeriod,Length,signalChoice[,C])
// 
// PARAMETERS:
//   Position ~ position of the sound source
//   MICS ~ microphones'' position
//   samplingPeriod ~ desired sampling period
//   Length ~ length of the signals
//   signal ~ signal function handle
//   C ~ sound speed
// 
// RETURN VALUE:
//   signals ~ generated discrete signals
//   times ~ sampling times
// 
// DESCRIPTION:
// 
//   The first syntax will generate as many signals as rows in MICS with 
//   the perceived signals from Position. The sampling period, the length 
//   and the signal should be specified. Sound speed is optional.
// 
//   The second syntax will generate the """"emitted"""" signal.
// 

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

// Input check
if %nargin<3 then
  error(strcat(["Usage: signals = GenerateDiscreteSignals(Position,MICS,samplingPeriod,Length,signal[,C])\n","Usage: signal = GenerateDiscreteSignals(samplingPeriod,Length,signal[,C])"]));
end;
if %nargin<5 then
  monoMode = %t;
  if %nargin<4 then
    C = 343.2;
  end;
else
  monoMode = %f;
  if %nargin<6 then
    C = 343.2;
  end;
end;

// Monomode?
if monoMode then
  // Remaping
  signal = samplingPeriod;
  Length = MICS;
  samplingPeriod = Position;
  // Times
  times = mtlb_imp(0,samplingPeriod,Length);
  // Computing
  signals = mtlb_e(signal,times);
else
  // Sampling times
  times = mtlb_imp(0,samplingPeriod,Length);
  // General variables
  // Number of microphones
  NMics = size(MICS,1);
  // Times of arrival
  TOA = sqrt(sum(mtlb_s(MICS,mtlb_repmat(Position,NMics,1)) .^2,2))/C;

  // Compute signals
  signals = zeros(NMics,length(times));
  for mic = 1:NMics
    signals(mic,:) = signal(times-TOA(mic));
  end;
end;

return;


endfunction
