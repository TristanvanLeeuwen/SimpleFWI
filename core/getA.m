function A = getA(f,m,h,n)
% Define 2D Helmholtz matrix with absorbing boundary conditions
%
% use:
%   A = getA(f,m,h,n)
%
% input:
%   f - frequency [Hz]
%   m - squared-slownes [s^2/km^2]
%   h - gridspacing in each direction d = [d1, d2];
%   n - number of gridpoints in each direction n = [n1, n2]
%
% output:
%   A - sparse matrix
%

%%
omega = 1e-3*2*pi*f;
N     = prod(n);

D1 = spdiags(ones(n(1),1)*[-1 1]/h(1),[0:1],n(1)-1,n(1));
D2 = spdiags(ones(n(2),1)*[-1 1]/h(2),[0:1],n(2)-1,n(2));
L  = -kron(speye(n(2)),D1'*D1) - kron(D2'*D2,speye(n(1)));

a = ones(n); a(:,[1 end]) = .5; a([1 end],:) = .5; a = a(:);

A = omega^2*diags(a.*m) + (2i*omega/h(1))*diags((1-a).*sqrt(m)) + L;