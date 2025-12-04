function [ica,Ae,cerr]=coshFpDeIca(X,epsilon,maxNumIter,debug)
% -------------------------------------------------------------------------------------------
% Fixed point algorithm (FastICA, Hyvarinen) using hyperbolic cosine (log(cosh())
% 
%  function [ica,Ae,cerr]=coshFpDeIca(X,epsilon,maxNumIter,debug)
%    X          : mixed signals (each row is a signal, must be centered)
%    epsilon    : convergence stopping criterion (default=0.0001)
%    maxNumIter : maximum number of iterations (default=200)
%    debug      : display each iteration of algorithm
%
%    ica   : separated signals 
%    Ae    : estimated mixing matrix (each row is a component)
%    cerr  : convergence error (1-> convergence failure)
%
% Simple version using as many variables as sources
% and deflationay orthogonalization
% 
% -------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% -------------------------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or modify it under the terms
% of the GNU General Public License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
%
% This program is distributed "as is" and "as available" in the hope that it will be useful,
% but WITHOUT ANY WARRANTY of any kind; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% -------------------------------------------------------------------------------------------

if(nargin<2), epsilon= 0.0001; end
if(nargin<3), maxNumIter=200; end
if(nargin<4), debug=0; end
  
[numVars, numSamples] = size(X);
cerr = 0;
% -------------------------- Whitening matrix estimation -----------------------

Rx=(X*X')/numSamples;     % covariance matrix
% Calculate the eigenvalues and eigenvectors of covariance matrix.
[E, D] = eig (Rx);
W = inv(sqrt (D)) * E';   % whitening matrix
% W = E*inv(sqrt (D))*E';   % whitening matrix  (W=sqrtm(Rx))
Z=W*X;     %  Whitened data


numOfIC=numVars;

a1=1;
fprintf('FastIca, cosh , deflationary \n');
fprintf('epsilon=%12.6f \n', epsilon);

W = zeros(numVars);
  
% The search for an independent component is repeated numOfIC times.
for ic = 1:numOfIC,
    
    fprintf('Component=%3d\n',ic);
    % Take a random initial vector
    w = rand(numVars, 1) - .5;
    w = w / norm(w);
 
%    wOld = w + epsilon;
    wOld = zeros(size(w));
    
    % Fixed-point iteration loop for a single IC.
    for iter = 1 : maxNumIter + 1
      
      % Test for termination condition. 
      % The algorithm converged if the direction of w and wOld is the same.
      absCos = abs(w' * wOld);
      if(1-absCos < epsilon), break, end
      if(debug), fprintf('it=%3d,%9.6f;  ',iter, 1-absCos);end
%       minNormDiff=min(norm(w - wOld), norm(w + wOld));
%       if minNormDiff < epsilon, break, end
%       fprintf('Iteration= %3d, %12.6f, eps=%12.6f \n',iter,  minNormDiff, epsilon);
      
      wOld=w;      
	    
      % Fixed point equation:
      % w <-- E[x g(w'x)] - E[g'(w'x)]w,   where g' is the derivative of g
      % chosing G(y)= 1/a log(cosh(ay)) we have: g=tanh(ay)) and  g'=a(1 -tanh^2(ay))
      % w <-- E[x tanh(a w'x)] - E[a (1-tanh^2(a w'x))] w
	  Tanh = tanh(a1 * Z' * w);
	  w = (Z * Tanh) / numSamples - a1 * sum(1 - Tanh .^ 2)'/ numSamples * w ;
     
      % Deflationary orthogonalization.
      % The current w vector is projected into the space orthogonal to the space  
      % spanned by the previuosly found W vectors. 
      w = w - W * W' * w;
      w = w / norm(w);
      
    end
    if(debug), fprintf('\n'); end
    fprintf('numIt=%3d,%9.6f\n',iter, 1-absCos);
    % Save the w vector
    W(:, ic) = w; 
    cerr = cerr | (1-absCos >= epsilon);
    if(cerr), fprintf('  ==> Convergence failure!\n');end
end
% if(cerr), Ae=[]; ica=[]; return; end

Ae = E*sqrt(D)*W;    % estimated mixing matrix

ica = W'*Z;    % separated components

end %== function ================================================================
%

