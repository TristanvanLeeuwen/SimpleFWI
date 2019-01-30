classdef opLInterp1D < opSpot
% 1D cubic Lagrange interpolation.
%
% use:
%   op = opLInterp1D(x1,x2)
%
% input:
%   x1 - input grid (i.e., [0:10:1000])
%   x2 - output grid (i.e., [0:1:1000])
%
% output:
%   op - SPOT operator: cubic Lagrange interpolation for input vectors of
%   same length as x1. Output has same lenght as x2.
%
% grids may also be non-equidistant and do not have to coincide
% completely. Output outside of the domain is set to zero.

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
        x1, x2;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opLInterp1D(x1,x2)
            op = op@opSpot('opLInterp1D', length(x2), length(x1));
            op.cflag     = 1;  
            op.linear    = 1; 
            op.children  = []; 
            op.sweepflag = true;
            op.x1        = x1;   
            op.x2        = x2;                      
            
        end %constructor
        
         
        function out = test(op)
            A  = getLA(op.x1,op.x2);
            % test1
            f1 = op.x1.^3;
            f2 = op.x2.^3;
            f2(op.x2<op.x1(1) | op.x2>op.x1(end)) = 0;
            e = norm(f2 - A*f1)/norm(f2);
            
            if e < 1e-10;
                %fprintf(1,'opLInterp1D: interpolation test succesfull: error = %g\n',e);
                out = 1;
            else
                %fprintf(2,'opLInterp1D: interpolation test failed    : error = %g\n',e);
                out  = 0;
            end
            
            % test2
            f1 = randn(size(f1));
            f2 = randn(size(f2));
            e = abs(1 - ((A*f1)'*f2) / (f1'*(A'*f2)));
            
            if e < 1e-10;
                %fprintf(1,'opLInterp1D: adjoint test succesfull: error = %g\n',e);
                out = out*1;
            else
                %fprintf(2,'opLInterp1D: adjoint test failed    : error = %g\n',e);
                out = out*0;
            end
            
        end
        
    end
    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = multiply(op,in,mode)
            A = getLA(op.x1,op.x2);
            if ~isempty(A)
                if mode==1
                    out = full(A*in);
                else
                    out = full(A'*in);
                end
            else
                fprintf(1,'Something wrong with grid\n');
                return;
            end
        end %multiply
       
    end %protected methods
    
end %classdef

function A = getLA(x1,x2)
    % find interior points
    ik = find(x2 >= x1(1) & x2 <= x1(end));
    % sizes
    n1=length(x1);n2=length(x2);nk=length(ik);
    % check
    if ~nk
        A = [];
        return;
    end
    
    % initialize stuff
    I = zeros(4*nk); J = I; S = I;
    a=1;b=2;c=3;d=4;
    l = 1;
    
    % loop
    for i = 1:nk
        k = ik(i);
        if x2(k)<x1(b)
            while (x2(k)<x1(b))&&(b-1>1)
                b=b-1;
            end
            a=b-1;c=b+1;d=c+1;
        elseif x2(k)>x1(c)
            while (x2(k)>x1(c))&&(c+1<n1)
                c=c+1;
            end
            a=c-2;b=c-1;d=c+1;
        end
        
        I(l)   = k; J(l)   = a; S(l)   = ((x2(k)-x1(b))*(x2(k)-x1(c))*(x2(k)-x1(d)))/((x1(a)-x1(b))*(x1(a)-x1(c))*(x1(a)-x1(d)));
        I(l+1) = k; J(l+1) = b; S(l+1) = ((x2(k)-x1(a))*(x2(k)-x1(c))*(x2(k)-x1(d)))/((x1(b)-x1(a))*(x1(b)-x1(c))*(x1(b)-x1(d)));
        I(l+2) = k; J(l+2) = c; S(l+2) = ((x2(k)-x1(b))*(x2(k)-x1(a))*(x2(k)-x1(d)))/((x1(c)-x1(b))*(x1(c)-x1(a))*(x1(c)-x1(d)));
        I(l+3) = k; J(l+3) = d; S(l+3) = ((x2(k)-x1(b))*(x2(k)-x1(c))*(x2(k)-x1(a)))/((x1(d)-x1(b))*(x1(d)-x1(c))*(x1(d)-x1(a)));
        l = l + 4;
    end
    % construct sparse matrix
    A = sparse(I(1:l-1),J(1:l-1),S(1:l-1),n2,n1); 
end


