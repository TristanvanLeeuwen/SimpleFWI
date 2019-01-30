function P = getP(h,n,zt,xt)
% Define sampling operator
%
% use:
%   P = getP(h,n,zt,xt)
%
% input:
%   h,n   - gridspacing and number of gridpoints
%   zt,xt - arrays defining sampling points (must coincide with grid)
%
% output
%   P     - spot operator


z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);

P = opLInterp2D(z(:), x(:), [zt(:), xt(:)]);