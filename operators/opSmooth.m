classdef opSmooth < opSpot
% Smoothing by convolving with [.25 .5 .25].
%
% use:
%   op = opSmooth(n,k)
%
% input:
%   n - length of input
%   k - number of times to convolve.
%
% output:
%   op - SPOT operator 

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: February, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        k;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opSmooth(n, k)

           op = op@opSpot('opSmooth',n,n);
           op.cflag     = 1;  
           op.linear    = 1;
           op.children  = []; 
           op.sweepflag = true;

           op.k         = k;      
       end %constructor
       
       function out = test(op)
       
           x  = randn(op.n,1);
           y  = mysmooth(x,op.k,1);
           xt = mysmooth(y,op.k,-1);
           
           e1 = abs((y'*y)/(x'*xt)-1);
           
           if~(e1<1e-10); fprintf(2,'opSmooth: adjoint test failed, error = %g\n',e1); end
           
           out = (e1<1e-10);
       end
       
    end
    
    
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function out = multiply(op,in,mode)
           if mode==1
               out = mysmooth(in,op.k,1);
           else
               out = mysmooth(in,op.k,-1);
           end
       end %multiply
    end %protected methods
    
end %classdef
   
    
function output=mysmooth(x,k,mode)
    n = size(x,1);
    f = [.25 .5 .25];

    S = spdiags(repmat(f,n,1),[-1:1],n,n); S(1,1)=S(1,1)+f(1); S(n,n)=S(n,n)+f(end);
    S = S^k;

    %convolve
    if mode==1
        output = full(S*x);
    else
        output = full(S'*x);
    end
end