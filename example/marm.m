%% setup

% read model, dx = 20 or 50
dx = 50;
v  = dlmread(['marm_' num2str(dx) '.dat']);

% initial model
v0 = @(zz,xx)v(1)+.7e-3*max(zz-350,0);

% set frequency, do not set larger than min(1e3*v(:))/(7.5*dx) or smaller
% than 0.5

f  = [1 2 3];

% receivers, xr = .1 - 10km, with 2*dx spacing, zr = 2*dx
xr = 100:2*dx:10000;
zr = 2*dx*ones(1,length(xr));

% sources, xr = .1 - 10km, with 4*dx spacing, zr = 2*dx
xs = 100:4*dx:10000;
zs = 2*dx*ones(1,length(xs));

% regularization parameter
alpha = 0;

%% observed data
% grid
n  = size(v);
h  = dx*[1 1];
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

% parameters
model.f = f;
model.n = n;
model.h = h;
model.zr = zr;
model.xr = xr;
model.zs = zs;
model.xs = xs;

% model
m = 1./v(:).^2;

% source
Q = eye(length(xs));

% data
D = F(m,Q,model);

%initial model
m0 = vec(1./v0(zz,xx).^2);

%% inversion
% misfit
fh = @(m)misfit(m,Q,D,alpha,model);

% Simple SD iteration
[mk,hist] = SDiter(fh,m0,1e-3,10,1e-6);

%% plot
vk = reshape(real(1./sqrt(mk)),n);

figure;
imagesc(x,z,v,[min(v(:)) max(v(:))]);title('ground truth');axis equal tight

figure;
imagesc(x,z,v0(zz,xx),[min(v(:)) max(v(:))]);title('initial');axis equal tight

figure;
imagesc(x,z,vk,[min(v(:)) max(v(:))]);title('final');axis equal tight

figure;
semilogy(hist(:,1),hist(:,2)/hist(1,2),hist(:,1),hist(:,3)/hist(1,3));title('convergence history');legend('f','|g|');