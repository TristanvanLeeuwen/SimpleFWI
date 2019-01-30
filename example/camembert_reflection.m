%% setup

% defined model
dx = 10;
v = 2 + 0.25*phantom([1, .5, .5, 0, 0, 0],1e3/dx + 1);

% initial model
v0 = @(z,x)2 + 0*z;

% set frequency, do not set larger than min(1e3*v(:))/(7.5*dx) or smaller
% than 0.5

f  = [5 10 15];

% receivers
xr = dx:2*dx:1e3-dx;
zr = 2*dx*ones(1,length(xr));

% sources
xs = dx:2*dx:1e3-dx;
zs = dx*ones(1,length(xr));

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

% data
D = F(m,model);

%initial model
m0 = vec(1./v0(zz,xx).^2);

%% inversion
% misfit
fh = @(m)misfit(m,D,alpha,model);

% Simple SD iteration
[mk,hist] = SDiter(fh,m0,1e-3,20,.1);

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