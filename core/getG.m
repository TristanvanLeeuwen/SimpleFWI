function G = getG(f,m,u,h,n)
% Define Jacobian matrix G(m,u) = d(A(m)*u)/dm
%
% use:
%   G = getG(f,m,u,d,n)
%
% input:
%   f - frequency [Hz]
%   m - squared-slownes [s^2/km^2]
%   u - wavefield
%   h - gridspacing in each direction d = [d1, d2];
%   n - number of gridpoints in each direction n = [n1, n2]
%
% output:
%   G - sparse matrix
%

%%
omega = 2*pi*f;
N     = prod(n);

% %
% w = [0 ones(1,n(1)-2) 0];
% if n(2)>1
%     w = w(:)*[0 ones(1,n(2)-2) 0];
% end
% w = w(:);
% 
% G = omega^2*spdiags(1e-6*w.*u,0,N,N) + 1i*omega*spdiags(1e-3*(1-w).*u./(2*sqrt(m)),0,N,N);

%%
omega = 1e-3*2*pi*f;
N     = prod(n);

a = ones(n); a(:,[1 end]) = .5; a([1 end],:) = .5; a = a(:);

G = omega^2*diags(a.*u) + (2i*omega/h(1))*diags(.5*(1-a).*u./sqrt(m));
