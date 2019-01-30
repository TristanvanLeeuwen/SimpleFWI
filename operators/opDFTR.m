classdef opDFTR < opSpot
% SPOT operator: DFT of real signal; output contains non-negative
% frequencies only. The adjoint is exact, at the expense of the inverse
% which differs roughly a factor two.
%
% use:
%   F = opDFTR(n)
%
% input:
%   n - length of input vector
%
% output:
%   F - FFT for real signals of length n, output will be of length floor(n/2) + 1
%

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
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opDFTR(n)
             
           m = floor(n/2) + 1;

           op = op@opSpot('opDFTR', m, n);
           op.cflag     = 1;  
           op.linear    = 1; 
           op.children  = [];    
           op.sweepflag = true;
           op.m         = m;    
           op.n         = n;                      
                                  
       end %constructor
       
       function out = test(op)
          x   = randn(op.n,10);
        
          xt1 = DFTR(op.m,op.n,x,1);
          xt2 = fft(x)/sqrt(op.n);
          xt3 = DFTR(op.m,op.n,xt1,-1);
          
          if mod(op.n,2)==1
                w = repmat([1; 2*ones(op.m-1,1)],1,10);
          else
                w = repmat([1; 2*ones(op.m-2,1); 1],1,10);
          end
          e1 = norm(xt1 - w.*xt2(1:op.m,:),'fro');
          %e2 = norm(x - xt3,'fro');

          x   = randn(op.n,1);
          y   = DFTR(op.m,op.n,x,1);
          e3 = abs((DFTR(op.m,op.n,x,1)'*y)/(x'*DFTR(op.m,op.n,y,-1)) - 1);
          
          if~(e1<1e-10); fprintf(2,'opDFTR: fft test failed, error = %g\n',e1); end
          if~(e3<1e-10); fprintf(2,'opDFTR: adjoint test failed, error = %g\n',e3);end
          
          out = (e1<1e-10)&(e3<1e-10);
       end
       
       
       
    end
    
    
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function out = multiply(op,in,mode)
           if mode==1
               out = DFTR(op.m,op.n,in,1);
           else
               out = DFTR(op.m,op.n,in,-1);
           end
           
       end %multiply
       
    end %protected methods
    
end %classdef


function output = DFTR(m,n,input,flag)

nc = size(input,2);
if flag == 1
    tmp = fft(input)/sqrt(n);
    if mod(n,2) == 1
        output = tmp(1:m,:) + [zeros(1,nc); conj(tmp(end:-1:m+1,:))];
    else
        output = tmp(1:m,:) + [zeros(1,nc); conj(tmp(end:-1:m+1,:)); zeros(1,nc)];
    end
else
    if mod(n,2) == 1
        tmp = [input ; conj(input(end:-1:2,:))];
    else
        tmp = [input ; conj(input(end-1:-1:2,:))];
    end
    output = ifft(tmp)*sqrt(n);
end
end
    
       