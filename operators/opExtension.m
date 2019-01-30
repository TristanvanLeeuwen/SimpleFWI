classdef opExtension < opSpot
% Extension operator. Pads input with constant values or zeros.
%
% use:
%   op = opExtension(n,nb,flag)
%
% input:
%   n    - length of input
%   nb   - number of points to add left nb(1) and right nb(2). If
%          length(nb)=1, use the same number of each side.
%   flag - 0: padd with zeros, 1: padd with boundary value.
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
    properties (SetAccess = private)
       nb;
       a;
    end % Properties


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opExtension(n,nb,a)
          if length(nb)==1
              nb = [nb nb];
          end
          nb = nb(1:2);
          m  = n + sum(nb);
          if nargin<3
              a = 1;
          end
          
          % Construct operator
          op = op@opSpot('Extension', m, n);
          op.nb = nb;
          op.sweepflag = 1;
          op.a = a;
       end % Constructor
       
       function out = test(op)
       
           x  = randn(op.n,1);
           y  = opExtension_intrnl(op.n,op.nb,op.a,x,1);
           xt = opExtension_intrnl(op.n,op.nb,op.a,y,-1);
           
           e1 = abs((y'*y)/(x'*xt)-1);
           
           if~(e1<1e-10); fprintf(2,'opExtension: adjoint test failed, error = %g\n',e1); end
           
           out = (e1<1e-10);
       end
       

    end % Methods
       
 
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
          y = opExtension_intrnl(op.n,op.nb,op.a,x,mode);
       end % Multiply          

    end % Methods
   
end % Classdef


%=======================================================================


function y = opExtension_intrnl(n,nb,a,x,mode)
    nx = size(x,2);
    if mode == 1
        if a
            y = [repmat(x(1,:),nb(1),1);x;repmat(x(end,:),nb(2),1)];
        else
            y = [zeros(nb(1),nx);x;zeros(nb(2),nx)];
        end
    else
       y = x(nb(1)+1:end-nb(2),:);
       if a
            y(1,:)   = y(1,:) + sum(x(1:nb(1),:),1);
            y(end,:) = y(end,:) + sum(x(end-nb(2)+1:end,:),1);
       end
    end
end

