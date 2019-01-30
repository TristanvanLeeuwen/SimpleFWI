classdef opSpline1D < opSpot
% 1D cubic spline evaluation.
%
% use:
%   op = opSpline1D(x1,x2,BC)
%
% input:
%   x1 - equidistant input grid
%   x2 - output grid, may be non-equidistant
%   BC - boundary conditions [left, right]. 0 = Dirichlet, 1 = Neumann

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
        x1,x2,BC;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opSpline1D(x1,x2,BC)
            op = op@opSpot('opSpline1D', length(x2), length(x1));
            op.cflag     = 1;  
            op.linear    = 1; 
            op.children  = []; 
            op.sweepflag = true;
  
            op.x1 = x1(:);    
            op.x2 = x2(:);
            op.BC = BC;
        end %constructor
        
        function out = test(op)
       
           x  = randn(op.n,1);
           y  = splneval1(op.x1,x,op.x2,[],1,op.BC);
           xt = splneval1(op.x1,[],op.x2,y,-1,op.BC);
           
           e1 = abs((y'*y)/(x'*xt)-1);
           
           if~(e1<1e-10); fprintf(2,'opSpline1D: adjoint test failed, error = %g\n',e1); end
           
           out = (e1<1e-10);
       end
    end
    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = multiply(op,in,mode)
            
            if mode==1
                out = splneval1(op.x1,in,op.x2,[],1,op.BC);
            else
                out = splneval1(op.x1,[],op.x2,in,-1,op.BC);
            end
            
            
        end %multiply
        
    end %protected methods
    
end %classdef

function output = splneval1(xi,yi,x,y,flag,bc)

h  = xi(2) - xi(1);
nh = 4;
xi = [xi(1) - h*[nh:-1:1]'; xi; xi(end) + h*[1:nh]'];
ni = length(xi);
nx = length(x);

A = [-1 3 -3 1; 3 -6 3 0; -3 0 3 0; 1 4 1 0]/6;

if flag>0
    m      = size(yi,2);
    output = zeros(nx,m);
    
    yi = [bc(1)*repmat(yi(1,:),nh,1); yi; bc(2)*repmat(yi(end,:),nh,1)];

    for k = 1:nx
        j = floor((x(k) - xi(1))/h) + 1;

        t = (x(k) - (xi(1) + (j-1)*h))/h;

        j = min(max(j,2),ni-2);

        output(k,:) = t.^[3:-1:0]*A*yi(j-1:j+2,:);
    end
else
    m      = size(y,2);
    output = zeros(ni,m);
    
    for k = 1:nx
        j = floor((x(k) - xi(1))/h) + 1;

        t = (x(k) - (xi(1) + (j-1)*h))/h;

        j = min(max(j,2),ni-2);

        output(j-1:j+2,:) = output(j-1:j+2,:) + (A'*t.^[3:-1:0]')*y(k,:);
    end
    
    output(nh+1,:)     = output(nh+1,:)   + bc(1)*sum(output(1:nh,:),1);
    output(end-nh,:)   = output(end-nh,:) + bc(2)*sum(output(end-nh+1:end,:),1);
    output             = output(nh+1:end-nh,:);
    
end

end


