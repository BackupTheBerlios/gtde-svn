function [J,GJ,HJ] = gTDECriterion(TDEs,PCCC,microphones,samplingPeriod)

// Output variables initialisation (not found in input variables)
J=[];
GJ=[];
HJ=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


//gTDECriterion implements the criterion to optimize for geometric TDE
// 
// USAGE: [J GJ HJ] = gTDECriterion(TDEs,PCCC,microphones,samplingPeriod)
// 
// PARAMETERS:
//   TDEs ~ set(s) of delays in which we want to evaluate the criterion
//   PCCC ~ the cross-correlation of the polynomial interpolation coefficients,
//     output of the PolynomialCrossCorrelation function.
//   microphones ~ positions of the microphones
//   samplingPeriod ~ samplingPeriod of the signal
// 
// RETURN VALUE:
//   J, GJ, HJ ~ criterion value, its gradient and hessian.
// 
// DESCRIPTION:
//     This function computes the determinant of the matrix of normalized
//     cross-correlation coefficients of the signals dealied accordingly to
//     TDEs.
// 
// REFERENCES:
//     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
//     estimation-based sound source localisation (gTDESSL). Research Report 
//     RR-7988, INRIA, June 2012.
// 
//   see also CCInterpolation, TDEis

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

//%% Input check
if %nargin<4 then
  error("Usage: [J GJ HJ] = gTDECriterion(TDEs, PCCC, microphones, samplingPeriod)");
end;

//%% General variables & check
// Number of microphones
NMics = size(microphones,1);
// TDEs dimension and number
[Dimension,NTDEs] = size(TDEs);

// Dimension vs NMics
if Dimension~=(NMics-1) then
  error("The dimesion of TDEs should be NMics-1 .");
end;

TDEs = TDEs*samplingPeriod;

//%% Auxiliar variables
// Vector r
r = zeros(Dimension,NTDEs);
// Vector rd
// % %     =========================================
// % %     =========================================
// % %     Be careful with the matrix rd. Each column does not represent the 
// % %     derivative of the vector r, but the derivative of the cross-correlation
// % %     function. Need to change the signs when using rd.
// % %     =========================================
// % %     =========================================
rd = zeros(Dimension,Dimension,NTDEs);
// Vector rdd
rdd = zeros(Dimension,Dimension,NTDEs);

// Fill the three vectors
for ss = 2:NMics
  // !! L.83: Unknown function TDEis not converted, original calling sequence used.
  index = TDEis(1,ss,NMics);
  // The cross-correlation function should be interpolated at
  // -t_(1,ss).entries not at t_(1,ss).entries
  // !! L.86: Unknown function CrossCorrelationInterpol not converted, original calling sequence used.
  
  [%v0,%v1,%v2] = CrossCorrelationInterpol(PCCC(1,ss).entries,-TDEs(index,:),samplingPeriod);
  r(index,:) = %v0;
  rd(index,index,:) = %v1;
  rdd(index,index,:) = %v2;
end;

// Matrix R
R = zeros(Dimension,Dimension,NTDEs);
// Matrix RD
RD = zeros(Dimension,Dimension,NTDEs);
// Matrix RDD
RDD = zeros(Dimension,Dimension,NTDEs);

// Fill the three matrices
for mic1 = 2:NMics
  // Autocorrelation
  // !! L.99: Unknown function CrossCorrelationInterpol not converted, original calling sequence used.
  R(mic1-1,mic1-1,:) = CrossCorrelationInterpol(PCCC(mic1,mic1).entries,0,samplingPeriod);
  for mic2 = mic1+1:NMics
    // Crosscorrelation
    // !! L.102: Unknown function TDEis not converted, original calling sequence used.
    // !! L.102: Unknown function TDEis not converted, original calling sequence used.
    TDEmic1mic2 = mtlb_s(TDEs(TDEis(1,mic2,NMics),:),TDEs(TDEis(1,mic1,NMics),:));
    // The cross-correlation function should be interpolated at
    // -t_(mic1,mic2).entries not at t_(mic1,mic2).entries
    // !! L.105: Unknown function CrossCorrelationInterpol not converted, original calling sequence used.
    
    [%v3,%v4,%v5] = CrossCorrelationInterpol(PCCC(mic1,mic2).entries,-TDEmic1mic2,samplingPeriod);
    R(mic1-1,mic2-1,:) = %v3;
    RD(mic1-1,mic2-1,:) = %v4;
    RDD(mic1-1,mic2-1,:) = %v5;
    // Symmetrize the matrices
    R(mic2-1,mic1-1,:) = R(mic1-1,mic2-1,:);
    RD(mic2-1,mic1-1,:) = RD(mic1-1,mic2-1,:);
    RDD(mic2-1,mic1-1,:) = RDD(mic1-1,mic2-1,:);
  end;
end;

//%% Compute the output

// 0. RTilde matrix
RTilde = zeros(NMics,NMics,NTDEs);
// !! L.117: Unknown function CrossCorrelationInterpol not converted, original calling sequence used.
RTilde(1,1,:) = CrossCorrelationInterpol(PCCC(1,1).entries,0,samplingPeriod);
RTilde(1,2:$,:) = r;
RTilde(2:$,1,:) = r;
RTilde(2:$,2:$,:) = R;

//%%% HERE WE NEED TO MAKE A DIFFERENCE BETWEEN NTDEs=1 and the rest
if NTDEs==1 then
  // 1. Criterion
  J = det(RTilde)/mtlb_prod(mtlb_diag(RTilde));

  // 2. If asked, compute the gradient 
  if %nargout>1 then
    // 2.1. Declare the gradient
    GJ = zeros(Dimension,1);
    // 2.2. Declare de auxiliar variable which stores the inverse or
    // RTilde times the first derivative of RTilde respecto to t1k
    RTildeID = zeros(NMics,NMics,Dimension);
    for tde = 1:Dimension
      // 2.3. First compute the derivative of the inverse of R
      mask = zeros(Dimension,Dimension);
      mask(tde,1:tde-1) = -1;
      mask(1:tde-1,tde) = -1;
      mask(tde,tde+1:$) = 1;
      mask(tde+1:$,tde) = 1;
      // 2.4. The derivative or RTilde
      RTildeD = zeros(NMics,NMics);
      RTildeD(1,2:$) = -rd(:,tde);
      RTildeD(2:$,1) = -rd(:,tde);
      RTildeD(2:$,2:$) = RD .*mask;
      // Storing
      RTildeID(:,:,tde) = mtlb_l(RTilde,RTildeD);
      // 2.5. Compute the mic-1''th derivative
      // !! L.149: Matlab function trace not yet converted, original calling sequence used.
      GJ = mtlb_i(GJ,tde,J*trace(RTildeID(:,:,tde)));
    end;
    GJ = GJ*samplingPeriod;
  end;
  // 3. If asked, compute the Hessian
  if %nargout>2 then
    // 3.1 Declare the Hessian
    HJ = zeros(Dimension,Dimension);
    // Fill the hessian
    for tde1 = 1:Dimension
      for tde2 = 1:Dimension
        // 3.1. The diagonal is different
        if tde1==tde2 then
          // 3.1.1. Compute the second derivative of R
          mask = zeros(Dimension,Dimension);
          mask(tde1,1:tde1-1) = 1;
          mask(1:tde1-1,tde1) = 1;
          mask(tde1,tde1+1:$) = 1;
          mask(tde1+1:$,tde1) = 1;
          auxRDD = RDD .*mask;
          // 3.1.2 Build the second derivative of RTilde
          RTildeDD = zeros(NMics,NMics);
          RTildeDD(1,2:$) = rdd(:,tde1);
          RTildeDD(2:$,1) = rdd(:,tde1);
          RTildeDD(2:$,2:$) = auxRDD;
          // 3.1.3. Compute the value
        
          // !! L.176: Matlab function trace not yet converted, original calling sequence used.
          // !! L.176: Matlab function trace not yet converted, original calling sequence used.
          HJ(tde1,tde2) = J*mtlb_a(trace(RTildeID(:,:,tde1)) .^2,trace(mtlb_a(-RTildeID(:,:,tde1)^2,mtlb_l(RTilde,RTildeDD))));
        else
          // 3.2. Not diagonal part
          // 3.2.1. Compute the second derivative of R
          mask = zeros(Dimension,Dimension);
          mask(tde1,tde2) = -1;
          mask(tde2,tde1) = -1;
          auxRDD = RDD .*mask;
          // 3.2.2. Compute the second derivative or RTilde
          RTildeDD = zeros(NMics,NMics);
          RTildeDD(2:$,2:$) = auxRDD;
          // 3.2.3. Add the remaining
        
          // !! L.189: Matlab function trace not yet converted, original calling sequence used.
          // !! L.189: Matlab function trace not yet converted, original calling sequence used.
          // !! L.189: Matlab function trace not yet converted, original calling sequence used.
          HJ(tde1,tde2) = J*mtlb_a(trace(RTildeID(:,:,tde2))*trace(RTildeID(:,:,tde1)),trace(mtlb_a(-RTildeID(:,:,tde2)*RTildeID(:,:,tde1),mtlb_l(RTilde,RTildeDD))));
        end;
      end;
    end;
    HJ = HJ*(samplingPeriod^2);
  end;
else
  //%%% Here NTDEs > 1!!!!!

  // 1. Criterion
  //criterionFun = @(M) det(M)/prod(diag(M));
  //J = squeeze(cellfun( criterionFun, mat2cell(RTilde,NMics,NMics,ones(NTDEs,1))))'';

  //         % Precompute the inverse matrices
  //         RTildeInverse = InverseMatrixArray(RTilde);

  // 2. If asked, compute the gradient 
  if %nargout>1 then
    // 2.1. Declare the gradient
    GJ = zeros(Dimension,NTDEs);
    // 2.2. Declare de auxiliar variable which stores the inverse or
    // RTilde times the first derivative of RTilde respecto to t1k
    RTildeID = zeros(NMics,NMics,Dimension,NTDEs);
    for tde = 1:Dimension
      // 2.3. First compute the derivative of the inverse of R
      mask = zeros(Dimension,Dimension,NTDEs);
      mask(tde,1:tde-1,:) = -1;
      mask(1:tde-1,tde,:) = -1;
      mask(tde,tde+1:$,:) = 1;
      mask(tde+1:$,tde,:) = 1;
      // 2.4. The derivative or RTilde
      RTildeD = zeros(NMics,NMics,NTDEs);
      RTildeD(1,2:$,:) = -rd(:,tde,:);
      RTildeD(2:$,1,:) = -rd(:,tde,:);
      RTildeD(2:$,2:$,:) = RD .*mask;
      for nt = 1:NTDEs
        // Storing
        RTildeID(:,:,tde,nt) = mtlb_l(RTilde(:,:,nt),RTildeD(:,:,nt));
        // 2.5. Compute the mic-1''th derivative
        // !! L.228: Matlab function squeeze not yet converted, original calling sequence used.
        // !! L.228: Matlab function trace not yet converted, original calling sequence used.
        GJ(tde,nt) = J(nt)*trace(squeeze(RTildeID(:,:,tde,nt)));
      end;
    end;
    GJ = GJ*samplingPeriod;
  end;

  // 3. If asked, compute the Hessian
  if %nargout>2 then
    // 3.1 Declare the Hessian
    HJ = zeros(Dimension,Dimension,NTDEs);
    // Fill the hessian
    for tde1 = 1:Dimension
      for tde2 = 1:Dimension
        // 3.1. The diagonal is different
        if tde1==tde2 then
          // 3.1.1. Compute the second derivative of R
          mask = zeros(Dimension,Dimension,NTDEs);
          mask(tde1,1:tde1-1,:) = 1;
          mask(1:tde1-1,tde1,:) = 1;
          mask(tde1,tde1+1:$,:) = 1;
          mask(tde1+1:$,tde1,:) = 1;
          auxRDD = RDD .*mask;
          // 3.1.2 Build the second derivative of RTilde
          RTildeDD = zeros(NMics,NMics,NTDEs);
          RTildeDD(1,2:$,:) = rdd(:,tde1,:);
          RTildeDD(2:$,1,:) = rdd(:,tde1,:);
          RTildeDD(2:$,2:$,:) = auxRDD;
          // 3.1.3. Compute the value
          for nt = 1:NTDEs
          
            // !! L.258: Matlab function trace not yet converted, original calling sequence used.
            // !! L.258: Matlab function squeeze not yet converted, original calling sequence used.
            // !! L.258: Matlab function squeeze not yet converted, original calling sequence used.
            // !! L.258: Matlab function trace not yet converted, original calling sequence used.
            HJ(tde1,tde2,nt) = J(nt)*mtlb_a(trace(RTildeID(:,:,tde1,nt)) .^2,trace(mtlb_a(-RTildeID(:,:,tde1,nt)^2,mtlb_l(squeeze(RTilde(:,:,nt)),squeeze(RTildeDD(:,:,nt))))));
          end;
        else
          // 3.2. Not diagonal part
          // 3.2.1. Compute the second derivative of R
          mask = zeros(Dimension,Dimension,NTDEs);
          mask(tde1,tde2,:) = -1;
          mask(tde2,tde1,:) = -1;
          auxRDD = RDD .*mask;
          // 3.2.2. Compute the second derivative or RTilde
          RTildeDD = zeros(NMics,NMics,NTDEs);
          RTildeDD(2:$,2:$,:) = auxRDD;
          // 3.2.3. Add the remaining
          for nt = 1:NTDEs
          
            // !! L.273: Matlab function trace not yet converted, original calling sequence used.
            // !! L.273: Matlab function trace not yet converted, original calling sequence used.
            // !! L.273: Matlab function trace not yet converted, original calling sequence used.
            HJ(tde1,tde2,nt) = J(nt)*mtlb_a(trace(RTildeID(:,:,tde2,nt))*trace(RTildeID(:,:,tde1,nt)),trace(mtlb_a(-RTildeID(:,:,tde2,nt)*RTildeID(:,:,tde1,nt),mtlb_l(RTilde(:,:,nt),RTildeDD(:,:,nt)))));
          end;
        end;
      end;
    end;
    HJ = HJ*(samplingPeriod^2);
  end;

end;

return;
endfunction
