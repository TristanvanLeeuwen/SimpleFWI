%% setup

% computational domain
dx = 20;
L = [1000, 1000];

% defined model
v = @(z,x) 2 + 0.25*((z-L(1)/2).^2 + (x-L(2)/2).^2 < (min(L)/4).^2);

% initial model
v0 = @(z,x)2 + 0*z;

% set frequency, do not set larger than min(1e3*v(:))/(7.5*dx) or smaller
% than 0.5

f = [2 3 4];
w = [.25,1,.25];

% receivers
nrec = 10;
zr = linspace(.1*L(1), .9*L(1), nrec);
xr = .1*L(2)*ones(1,nrec);

% sources
nsrc = 10;
zs = linspace(.1*L(1), .9*L(1), nsrc);
xs = .9*L(2)*ones(1,nsrc);

% regularization parameter
alpha = 0;

%% observed data

% grid
n  = L/dx + 1;
h  = dx*[1 1];
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

% parameters
model.f = f;
model.w = w;
model.n = n;
model.h = h;
model.zr = zr;
model.xr = xr;
model.zs = zs;
model.xs = xs;

% model
m = 1./v(zz(:),xx(:)).^2;

% data
Q = eye(length(xs));
D = F(m,Q,model);

%initial model
m0 = vec(1./v0(zz,xx).^2);

%% inversion
% misfit
fh = @(m)misfit(m,Q,D,alpha,model);

% Simple SD iteration
[mk,hist] = SDiter(fh,m0,1e-3,20,1e-6);

%% plot
vk = reshape(real(1./sqrt(mk)),n);
vt = reshape(real(1./sqrt(m)),n);

figure;
imagesc(x,z,vt,[min(vt(:)) max(vt(:))]);title('ground truth');axis equal tight

figure;
imagesc(x,z,v0(zz,xx),[min(vt(:)) max(vt(:))]);title('initial');axis equal tight

figure;
imagesc(x,z,vk,[min(vt(:)) max(vt(:))]);title('final');axis equal tight

figure;
semilogy(hist(:,1),hist(:,2)/hist(1,2),hist(:,1),hist(:,3)/hist(1,3));title('convergence history');legend('f','|g|');