classdef opLInterp2D < opSpot
% 2D linear Lagrange interpolation.
%
% use:
%   op = opLInterp2D(x1,y1,X2)
% 
% input:
%   x1 (fast dim),y1 (slow dim) - rectangular input grid (i.e., [0:10:1000],[0:10:1000])
%   X2 = [x2 y2]                - list of evaluation points in N x 2 matrix
%
% input grid may also be non-equidistant. Output outside of the domain is set to zero.
% It's faster to use a Kronecker of opLInterp1D when evaluating on a
% rectangular grid.

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
        %put any variables that you need to store here
        % for instance, any variables you pass to the function handle
        % should be stored here instead, and then they will be available in
        % the multiply function
        % m, n, linear, cflag, children - are already provided by
        % superclass no need to redefine them here.
        x1, y1, X2;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opLInterp2D(x1,y1,X2)
  
            op = op@opSpot('opLInterp2D', size(X2,1), length(x1)*length(y1));
            op.cflag     = 1;  
            op.linear    = 1; 
            op.children  = [];
            
            op.x1 = x1;    
            op.y1 = y1;
            op.X2  = X2;
        end %constructor
        
        function out = test(op)
            A  = getLA2(op.x1,op.y1,op.X2);
            
            % test1
            [xx,yy] = ndgrid(op.x1,op.y1);
            f1 = xx(:) + 2*yy(:);
            f2 = op.X2(:,1) + 2*op.X2(:,2);
            
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
            A  = getLA2(op.x1,op.y1,op.X2);
            if ~isempty(A)
                if mode==1
                    out = A*in;
                else
                    out = A'*in;
                end
            else
                fprintf(1,'Something wrong with grid\n');
                return;
            end
            
        end %multiply
        
    end %protected methods
    
end %classdef

function A = getLA2(x1,y1,X2)
    % interior points
    ik = find(X2(:,1) >= x1(1) & X2(:,1) <= x1(end) & X2(:,2) >=y1(1) & X2(:,2) <= y1(end));
    % sizes
    nx1 = length(x1); ny1 = length(y1);
    n1  = nx1*ny1;
    n2  = size(X2,1);
    nk  = length(ik);
    % check
    if ~nk
        A = [];
        return;
    end
    % index sets for sparse matrix
    I = zeros(4*nk); J = I; S = J;
    % loop
    l = 1;
    for i = 1:nk
        k = ik(i);
        ix = min(find(x1<=X2(k,1), 1, 'last' ), nx1 - 1);
        iy = min(find(y1<=X2(k,2), 1, 'last' ), ny1 - 1);
        a = ix     + nx1*(iy - 1);
        b = ix     + nx1*(iy);
        c = ix + 1 + nx1*(iy - 1);
        d = ix + 1 + nx1*(iy);
        
        I(l)   = k; J(l)   = a; S(l)   = (X2(k,1) - x1(ix+1))*(X2(k,2) - y1(iy+1))/((x1(ix) - x1(ix+1))*(y1(iy) - y1(iy+1)));
        I(l+1) = k; J(l+1) = b; S(l+1) = (X2(k,1) - x1(ix+1))*(X2(k,2) - y1(iy))  /((x1(ix) - x1(ix+1))*(y1(iy+1) - y1(iy)));
        I(l+2) = k; J(l+2) = c; S(l+2) = (X2(k,1) - x1(ix))  *(X2(k,2) - y1(iy+1))/((x1(ix+1) - x1(ix))*(y1(iy) - y1(iy+1)));
        I(l+3) = k; J(l+3) = d; S(l+3) = (X2(k,1) - x1(ix))  *(X2(k,2) - y1(iy))  /((x1(ix+1) - x1(ix))*(y1(iy+1) - y1(iy)));
        l = l + 4;
    end
    % construct matrix
    A = sparse(I(1:l-1),J(1:l-1),S(1:l-1),n2,n1);
end





