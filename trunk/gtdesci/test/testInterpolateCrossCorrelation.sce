
// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


// test the interpolation of the cross-correlation function

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

xdel(winsid());clear;clc

getd('../src/');

// Scale parameters
a = 0.1;
c = a/sqrt(2);
// Microphones in a regular tetrahedron
MICS = [a,0,-c;-a,0,-c;0,a,c;0,-a,c];

// Position
position = [1.34, 2.3, -2.4];

// Sampling frequency
samplingPeriod = 1/48000;

// Length
Length = 0.01;

// Generate the signal
signals = GenerateDiscreteSignals(position,MICS,samplingPeriod,Length,testSignalSCILAB);

// Plotting period
plottingPeriod = samplingPeriod/50;
h = 0.9*plottingPeriod;
// Plotting times
// !! L.70: Unknown function TDEmax not converted, original calling sequence used.
maxTDE = TDEmax(MICS);
maxTDE = mtlb_e(maxTDE,1);
maxLAG = ceil(maxTDE/samplingPeriod);
maxTDE = maxLAG*samplingPeriod;
plottingTimes = mtlb_imp(-maxTDE,plottingPeriod,maxTDE);
// sampling times
samplingTimes = 0:samplingPeriod:Length;
// !! L.77: Unknown function InterpolateCrossCorrelat not converted, original calling sequence used.
[myCorr,myDerivative,my2Derivative] = InterpolateXCorrelation(signals(1,:),signals(2,:),plottingTimes,samplingPeriod,maxTDE);
//     myCorr = myCorr / max(myCorr);
//     myDerivative = myDerivative / max(myDerivative);
//     my2Derivative = my2Derivative / max(my2Derivative);

// Discrete correlation
//     [discreteCorr lags] = xcorr(signals(1,:),signals(2,:),maxLAG,''unbiased'');
//     discreteCorr = discreteCorr / max(discreteCorr);

t1 = norm(Positions(positionIndex,:)-MICS(1,:),2)/343.2;
t2 = norm(Positions(positionIndex,:)-MICS(2,:),2)/343.2;
// L.88: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf("TDE: %1.8e\n",-(t2-t1));

aCC = analyticCrossCorrelation(F,1/L,Length,t1,t2,plottingTimes);
//     aCC = aCC / max(aCC);

// First plot: signals
subplot(2,2,1);
plot(samplingTimes,signals(1,:));
set(gca(),"auto_clear","off")
plot(samplingTimes,signals(2,:),"r");
set(gca(),"auto_clear","on")
// !! L.99: Matlab function legend not yet converted, original calling sequence used.
legend("signal1","signal2");

// Second plot: cross-correlation
subplot(2,2,2);
//     plot( -maxTDE:samplingPeriod:maxTDE, discreteCorr, ''ob'');
set(gca(),"auto_clear","off")
plot(plottingTimes,myCorr,"gx");
plot(plottingTimes,aCC,"k");
set(gca(),"auto_clear","on")
// !! L.108: Matlab function legend not yet converted, original calling sequence used.
legend("Mycorr","Analytic autocorrelation");

// Thir plot: xcorr''s derivative
subplot(2,2,3);
//     dDerivative = discreteCorr(2:end)-discreteCorr(1:end-1);
//     dDerivative = dDerivative / max(dDerivative);
//     plot( - maxTDE:samplingPeriod:(maxTDE-samplingPeriod), dDerivative,''ob'');
set(gca(),"auto_clear","off")
plot(plottingTimes,myDerivative,"gx");
set(gca(),"auto_clear","on")
// !! L.118: Matlab function legend not yet converted, original calling sequence used.
legend("My derivative");

// Fourth: 2nd derivative
subplot(2,2,4);
//     ddDerivative = dDerivative(2:end)-dDerivative(1:end-1);
//     ddDerivative = ddDerivative / max(ddDerivative);
//     plot( - maxTDE:samplingPeriod:(maxTDE-2*samplingPeriod), ddDerivative,''ob'');
set(gca(),"auto_clear","off")
plot(plottingTimes,my2Derivative,"gx");
set(gca(),"auto_clear","on")
// !! L.128: Matlab function legend not yet converted, original calling sequence used.
legend("My 2nd derivative");
return;S

function [acc] = analyticCrossCorrelation(F,lambda,%length,t1,t2,tau)
    // Output variables initialisation (not found in input variables)
    acc=[];
    
    // Display mode
    mode(0);
    
    // Display warning for floating point exception
    ieee(1);
    
    
    M = mtlb_min(mtlb_a(%length,t1),mtlb_a(mtlb_a(%length,t2),tau));
    m = mtlb_max(t1,mtlb_a(t2,tau));
    
    
    acc = ((lambda/4)*exp(mtlb_a(mtlb_a(t1,t2),tau)/lambda)) .*mtlb_s(mtlb_s(exp(-(2*M)/lambda) .*mtlb_s(cos(F*mtlb_s(mtlb_s(mtlb_s(2*M,t1),tau),t2)),(F*lambda)*sin(F*mtlb_s(mtlb_s(mtlb_s(2*M,t1),tau),t2))),exp(-(2*m)/lambda) .*mtlb_s(cos(F*mtlb_s(mtlb_s(mtlb_s(2*m,t1),t2),tau)),(F*lambda)*sin(F*mtlb_s(mtlb_s(mtlb_s(2*m,t1),t2),tau))))/mtlb_a(1,(lambda*F) .^2),cos(F*mtlb_s(mtlb_a(t2,tau),t1)) .*mtlb_s(exp(-(2*M)/lambda),exp(-(2*m)/lambda)));
    
    return;
endfunction
