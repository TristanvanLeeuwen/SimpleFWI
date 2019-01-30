function [D,J] = F(m,model)
% Forward operator
%
%   D = P^TA^{-1}(m)Q
%
% where P, Q encode the receiver and source locations and L is the first-order FD matrix
%
% use:
%   [D,J] = F(m,model);
%
% input:
%   m - squared-slownes [s^2/km^2]
%   model.h - gridspacing in each direction d = [d1, d2];
%   model.n - number of gridpoints in each direction n = [n1, n2]
%   model.f - frequency [Hz].
%   model.{zr,xr} - {z,x} locations of receivers [m] (must coincide with gridpoints)
%   model.{zs,xs} - {z,x} locations of sources [m] (must coincide with gridpoints)
%
%
% output:
%   D - data matrix
%   J - Jacobian as Spot operator

    % size
    nr = length(model.zr);
    ns = length(model.zs);
    nf = length(model.f);
    nx = prod(model.n);
    
    % generate matrices
    P = getP(model.h,model.n,model.zr,model.xr);
    Q = getP(model.h,model.n,model.zs,model.xs);

    % solve
    U = zeros(nx,ns,nf);
    D = zeros(nr,ns,nf);
    for k = 1:nf
        Ak = getA(model.f(k),m,model.h,model.n);
        U(:,:,k) = Ak\Q;
        D(:,:,k) = P'*U(:,:,k);
    end
    D = D(:);
    % Jacobian
    J = opFunction(nr*ns*nf, nx, @(x,flag)Jmv(x,m,U,model,flag));
end


function y = Jmv(x,m,U,model,flag)
    % size
    nr = length(model.zr);
    ns = length(model.zs);
    nf = length(model.f);
    nx = prod(model.n);
    
    %% get matrices
    P = getP(model.h,model.n,model.zr,model.xr);
 
    %% compute mat-vec
    if flag == 1
        y = zeros(nr,ns,nf);
        for k = 1:nf
            Rk = zeros(nx,ns);
            Ak = getA(model.f(k),m,model.h,model.n);
            Gk = @(u)getG(model.f(k),m,u,model.h,model.n);
            for l = 1:ns
               Rk(:,l) = -Gk(U(:,l,k))*x;
            end
            y(:,:,k) = P'*(Ak\Rk);
        end
        y = y(:);
    else
        y = zeros(nx,1);
        x = reshape(x,[nr,ns,nf]);
        for k = 1:nf
            Ak = getA(model.f(k),m,model.h,model.n);
            Gk = @(u)getG(model.f(k),m,u,model.h,model.n);
            Rk = Ak'\(P*x(:,:,k));
            for l = 1:size(U,2)
                y = y - Gk(U(:,l,k))'*Rk(:,l);
            end
        end
    end

end
