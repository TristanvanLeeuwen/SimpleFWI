classdef opSpline2D < opSpot
% 2D cubic spline evaluation.
%
% use:
%   op = opSpline2D(x1,y1,x2,y2,BC)
%
% input:
%   x1 (fast dim),y1 (slow dim) - equidistant input grid
%   x2 (fast dim),y2 (slow dim) - output grid
%   BC - boundary conditions [top, bottom, left, right]. 0 = Dirichlet, 1 = Neumann.
%
% in principle this gives the same result as 
% opKron(opSpline1D(y1,y2,[left right]),opSpline1D(x1,x2,[top bottom]));

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
        x1,y1,x2,y2,BC;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opSpline2D(x1,y1,x2,y2,BC)
            op = op@opSpot('opSpline2D', length(x2)*length(y2), length(x1)*length(y1));
            op.cflag     = 1;  
            op.linear    = 1; 
            op.children  = [];
         
            op.x1 = x1;   
            op.x2 = x2;
            op.y1 = y1;
            op.y2 = y2;
            op.BC = BC;
        end %constructor
        
    end
    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = multiply(op,in,mode)
            
            if mode==1
                out = splneval2(op.x1,op.y1,in,op.x2,op.y2,[],1,op.BC);
            else
                out = splneval2(op.x1,op.y1,[],op.x2,op.y2,in,-1,op.BC);
            end
            
            
        end %multiply
        
    end %protected methods
    
end %classdef

function output = splneval2(xi,yi,si,x,y,s,flag,bc)
% 2D Cubic spline evaluation
% use:		output = lagintp(x, y, s, x1, y1, s1, flag,{bc})
% input:	x,y - grid
%           s - corresponding data, for forward mode
%           x1,y1 - grid
%           s1 - corresponding data, for adjoint mode
%           flag - 1:forward, 2:adjoint
%           bc - boundary conditions [top bot left right], 0=dirichlet,1=neumann.
%                default: bc=[1 1 1 1]
%output: 	forward: s1, adjoint: s
%note:		all vectors are assumed column vectors
%
if (size(xi)==size(x))&(size(yi)==size(y))
    if ~norm(xi-x)&~norm(yi-y)
        if flag>0
            output = si(:);
        else
            output = s(:);
        end
        return;
    end
end
if (nargin<8)|(length(bc)<4)
    bc = [1 1 1 1];
end

nx = length(x); dx = x(2)-x(1);
ny = length(y); dy = y(2)-y(1);

nxi = length(xi); dxi = xi(2)-xi(1);
nyi = length(yi); dyi = yi(2)-yi(1);

%extend grid
xi = [xi(1)-2*dxi:dxi:xi(end)+2*dxi]';
yi = [yi(1)-2*dyi:dyi:yi(end)+2*dyi]';
nxi = length(xi); dxi = xi(2)-xi(1);
nyi = length(yi); dyi = yi(2)-yi(1);
%
if flag==1
    
    si = reshape(si,nxi-4,nyi-4);
    output = zeros(nx,ny);
    
    %bc's
    %top,bottom
    si = [repmat(bc(1)*si(1,:),2,1);si;repmat(bc(2)*si(end,:),2,1)];
    %left,right
    si = [repmat(bc(3)*si(:,1),1,2) si repmat(bc(4)*si(:,end),1,2)];
    
    for k = 1:nx
        for l = 1:ny
            xh = (x(k)-xi(1))/dxi;
            ih = floor(xh); xr = xh-ih; ih = ih+1;
            im = max(1,min(nxi,ih-1));i0 = max(1,min(nxi,ih));
            i1 = max(1,min(nxi,ih+1));i2 = max(1,min(nxi,ih+2));
            %
            yh = (y(l)-yi(1))/dyi;
            jh = floor(yh); yr = yh-jh; jh = jh+1;
            jm = max(1,min(nyi,jh-1));j0 = max(1,min(nyi,jh));
            j1 = max(1,min(nyi,jh+1));j2 = max(1,min(nyi,jh+2));
            %
            cimjm = si(im,jm); ci0jm = si(i0,jm); ci1jm = si(i1,jm); ci2jm = si(i2,jm);
            cimj0 = si(im,j0); ci0j0 = si(i0,j0); ci1j0 = si(i1,j0); ci2j0 = si(i2,j0);
            cimj1 = si(im,j1); ci0j1 = si(i0,j1); ci1j1 = si(i1,j1); ci2j1 = si(i2,j1);
            cimj2 = si(im,j2); ci0j2 = si(i0,j2); ci1j2 = si(i1,j2); ci2j2 = si(i2,j2);
            %
            sh0 = (1./36.)*( 16*ci0j0+ 4*(cimj0+ci0jm+ci0j1+ci1j0)+ (cimjm+ci1jm+cimj1+ci1j1) );
            sh1 = (1./12.)*(4*(ci1j0-cimj0) + ((ci1j1-cimjm)+(ci1jm-cimj1) ));
            sh2 = (1./12.)*(4*(ci0j1-ci0jm) + ((ci1j1-cimjm)-(ci1jm-cimj1) ));
            sh3 = (1./12.)*( (ci1j1+ci1jm+cimj1+cimjm)+ 4*(ci1j0+cimj0)-2*(ci0jm+ci0j1)-8*ci0j0);
            sh4 = 0.25*( (ci1j1-ci1jm) + (cimjm-cimj1));
            sh5 = (1./12.)*( (ci1j1+ci1jm+cimj1+cimjm)+ 4*(ci0j1+ci0jm)-2*(ci1j0+cimj0)-8*ci0j0);
            sh6 = (1./36.)*(12*(ci0j0-ci1j0)+4*(ci2j0-cimj0)+3*((ci0j1-ci1j1)+(ci0jm-ci1jm))+((ci2j1-cimj1)+(ci2jm-cimjm)));
            sh7 = 0.25*(2*(ci0jm-ci0j1)+(ci1j1-ci1jm)+(cimj1-cimjm));
            sh8 = 0.25*(2*(cimj0-ci1j0)+(ci1j1-cimj1)+(ci1jm-cimjm));
            sh9 = (1./36.)*(12*(ci0j0-ci0j1)+4*(ci0j2-ci0jm)+3*((ci1j0-ci1j1)+(cimj0-cimj1))+((ci1j2-ci1jm)+(cimj2-cimjm)) );
            sh10 = (1./12.)*(3*((ci0j1-ci0jm)+(ci1jm-ci1j1))+(ci2j1-ci2jm)+(cimjm-cimj1));
            sh11 = ci0j0-0.5*(cimj0+ci0jm+ci0j1+ci1j0)+0.25*(cimjm+cimj1+ci1jm+ci1j1);
            sh12 = (1./12.)*(3*((ci1j0-ci1j1)+(cimj1-cimj0) )+(cimjm-ci1jm)+(ci1j2-cimj2) );
            sh13 = (1./12.)*(6*(ci1j0-ci0j0)+3*((ci0j1-ci1jm)+(ci0jm-ci1j1))+2*(cimj0-ci2j0)+(ci2jm-cimjm)+(ci2j1-cimj1) );
            sh14 = (1./12.)*(6*(ci0j1-ci0j0)+3*((cimj0-cimj1)+(ci1j0-ci1j1))+2*(ci0jm-ci0j2)+(ci1j2-cimjm)+(cimj2-ci1jm) );
            sh15 = (1./36.)*(9*((ci0j0-ci0j1)+(ci1j1-ci1j0))+3*((ci2j0-ci2j1)+(cimj1-cimj0)+(ci0j2-ci1j2)+(ci1jm-ci0jm))+((cimjm-ci2jm)+ (ci2j2-cimj2)) );
            %
            v0   = sh0 + xr*(sh1  + xr*(sh3  + xr*sh6 ));
            v1   = sh2 + xr*(sh4  + xr*(sh7  + xr*sh10));
            v2   = sh5 + xr*(sh8  + xr*(sh11 + xr*sh13));
            v3   = sh9 + xr*(sh12 + xr*(sh14 + xr*sh15));
            %
            output(k,l) = v0 + yr*(v1 + yr*(v2 + yr*v3));
        end
    end
elseif flag==-1
    output = zeros(nxi,nyi);
    s = reshape(s,nx,ny);
    for k = 1:nx
        for l = 1:ny
            xh = (x(k)-xi(1))/dxi;
            ih = floor(xh); xr = xh-ih; ih = ih+1;
            im = max(1,min(nxi,ih-1));i0 = max(1,min(nxi,ih));
            i1 = max(1,min(nxi,ih+1));i2 = max(1,min(nxi,ih+2));
            %
            yh = (y(l)-yi(1))/dyi;
            jh = floor(yh); yr = yh-jh; jh = jh+1;
            jm = max(1,min(nyi,jh-1));j0 = max(1,min(nyi,jh));
            j1 = max(1,min(nyi,jh+1));j2 = max(1,min(nyi,jh+2));
            
            ci2j2 = (1/36)*xr^3*yr^3;
            ci1j0 = ((1/9 + xr/3 + xr^2/3 - xr^3/3 - yr^2/6 - 0.5*xr*yr^2 - 0.5*xr^2*yr^2 + (xr^3*yr^2)/2 + yr^3/12 + (xr*yr^3)/4 + (xr^2*yr^3)/4 - (xr^3*yr^3)/4));
            ci0j1 = ((1/9 - xr^2/6 + xr^3/12 + yr/3 - 0.5*xr^2*yr + (xr^3*yr)/4 + yr^2/3 - 0.5*xr^2*yr^2 + (xr^3*yr^2)/4 - yr^3/3 + (xr^2*yr^3)/2 - (xr^3*yr^3)/4));
            ci2j1 = ((xr^3/36 + (xr^3*yr)/12 + (xr^3*yr^2)/12 - (xr^3*yr^3)/12));
            ci1j2 = ((yr^3/36 + (xr*yr^3)/12 + (xr^2*yr^3)/12 - (xr^3*yr^3)/12));
            ci0jm = ((1/9 - xr^2/6 + xr^3/12 - yr/3 + 0.5*xr^2*yr - (xr^3*yr)/4 + yr^2/3 - 0.5*xr^2*yr^2 + (xr^3*yr^2)/4 - yr^3/9 + (xr^2*yr^3)/6 - (xr^3*yr^3)/12));
            cimj0 = ((1/9 - xr/3 + xr^2/3 - xr^3/9 - yr^2/6 + 0.5*xr*yr^2 - 0.5*xr^2*yr^2 + (xr^3*yr^2)/6 + yr^3/12 - (xr*yr^3)/4 + (xr^2*yr^3)/4 - (xr^3*yr^3)/12));
            ci2jm = ((xr^3/36 - (xr^3*yr)/12 + (xr^3*yr^2)/12 - (xr^3*yr^3)/36));
            cimj2 = ((yr^3/36 - (xr*yr^3)/12 + (xr^2*yr^3)/12 - (xr^3*yr^3)/36));
            cimjm = ((1/36 - xr/12 + xr^2/12 - xr^3/36 - yr/12 + 0.25*xr*yr - 0.25*xr^2*yr + (xr^3*yr)/12 + yr^2/12 - 0.25*xr*yr^2 + 0.25*xr^2*yr^2 - (xr^3*yr^2)/12 - yr^3/36 + (xr*yr^3)/12 - (xr^2*yr^3)/12 + (xr^3*yr^3)/36));
            ci2j0 = ((xr^3/9 - (xr^3*yr^2)/6 + (xr^3*yr^3)/12));
            cimj1 = ((1/36 - xr/12 + xr^2/12 - xr^3/36 + yr/12 - 0.25*xr*yr + 0.25*xr^2*yr - (xr^3*yr)/12 + yr^2/12 - 0.25*xr*yr^2 + 0.25*xr^2*yr^2 - (xr^3*yr^2)/12 - yr^3/12 + (xr*yr^3)/4 - (xr^2*yr^3)/4 + (xr^3*yr^3)/12));
            ci0j2 = ((yr^3/9 - (xr^2*yr^3)/6 + (xr^3*yr^3)/12));
            ci1jm = ((1/36 + xr/12 + xr^2/12 - xr^3/12 - yr/12 - 0.25*xr*yr - 0.25*xr^2*yr + (xr^3*yr)/4 + yr^2/12 + 0.25*xr*yr^2 + 0.25*xr^2*yr^2 - (xr^3*yr^2)/4 - yr^3/36 - (xr*yr^3)/12 - (xr^2*yr^3)/12 + (xr^3*yr^3)/12));
            ci0j0 = ((4/9 - (2*xr^2)/3 + xr^3/3 - (2*yr^2)/3 + xr^2*yr^2 - (xr^3*yr^2)/2 +  yr^3/3 - (xr^2*yr^3)/2 + (xr^3*yr^3)/4));
            ci1j1 = ((1/36 + xr/12 + xr^2/12 - xr^3/12 + yr/12 + 0.25*xr*yr + 0.25*xr^2*yr - (xr^3*yr)/4 + yr^2/12 + 0.25*xr*yr^2 + 0.25*xr^2*yr^2 - (xr^3*yr^2)/4 - yr^3/12 - (xr*yr^3)/4 - (xr^2*yr^3)/4 + (xr^3*yr^3)/4));
            
            
            output(im,jm) = output(im,jm)+cimjm*s(k,l);
            output(im,j0) = output(im,j0)+cimj0*s(k,l);
            output(im,j1) = output(im,j1)+cimj1*s(k,l);
            output(im,j2) = output(im,j2)+cimj2*s(k,l);
            output(i0,jm) = output(i0,jm)+ci0jm*s(k,l);
            output(i0,j0) = output(i0,j0)+ci0j0*s(k,l);
            output(i0,j1) = output(i0,j1)+ci0j1*s(k,l);
            output(i0,j2) = output(i0,j2)+ci0j2*s(k,l);
            output(i1,jm) = output(i1,jm)+ci1jm*s(k,l);
            output(i1,j0) = output(i1,j0)+ci1j0*s(k,l);
            output(i1,j1) = output(i1,j1)+ci1j1*s(k,l);
            output(i1,j2) = output(i1,j2)+ci1j2*s(k,l);
            output(i2,jm) = output(i2,jm)+ci2jm*s(k,l);
            output(i2,j0) = output(i2,j0)+ci2j0*s(k,l);
            output(i2,j1) = output(i2,j1)+ci2j1*s(k,l);
            output(i2,j2) = output(i2,j2)+ci2j2*s(k,l);
            
        end
    end
    
    % left,right
    output(:,3) = output(:,3) + bc(3)*sum(output(:,1:2),2);
    output(:,end-2) = output(:,end-2) + bc(4)*sum(output(:,end-1:end),2);
    % top,bottom
    output(3,:) = output(3,:) + bc(1)*sum(output(1:2,:),1);
    output(end-2,:) = output(end-2,:) + bc(2)*sum(output(end-1:end,:),1);
    
    output = output(3:end-2,3:end-2);
end
output = output(:);
%%EOF
end