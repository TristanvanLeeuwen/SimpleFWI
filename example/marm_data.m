%% setup

% read model, dx = 20 or 50
dx = 20;
v  = dlmread(['marm_' num2str(dx) '.dat']);

% initial model
v0 = @(zz,xx)v(1)+.7e-3*max(zz-350,0);

% frequency
T = 10; % 10 seconds of data
dt = 0.001; % sampling
f = [0:1/T:30];

% wavelet
f0 = 10;
t0 = 0.1;
w = f.^2.*exp(-(f/f0).^2).*exp(1i*2*pi*f*t0);

% receivers, xr = .1 - 10km, with 2*dx spacing, zr = 2*dx
xr = 100:2*dx:10000;
zr = 2*dx*ones(1,length(xr));

% sources, xr = .1 - 10km, with 4*dx spacing, zr = 2*dx
xs = 100:4*dx:10000;
zs = 2*dx*ones(1,length(xs));

% regularization parameter
alpha = 0.1;

%% observed data
% grid
n  = size(v);
h  = dx*[1 1];
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

% parameters
model.n = n;
model.h = h;
model.zr = zr;
model.xr = xr;
model.zs = zs;
model.xs = xs;

% model
m = 1./v(:).^2;

%% data
D = zeros(length(xr),length(xs),length(f));
for k = 1:length(f)
    model.f = f(k);
    D(:,:,k) = w(k)*full(F(m,model));
end
D = permute(D,[3,1,2]);

%% Fourier transformatie
t = (0:dt:T)';
Dt = zeros(length(t),length(xr),length(xs));
for k = 1:length(f)
    %Dt = Dt + exp(-1i*2*pi*f(k)*t)*reshape(D(:,:,k),1,[]) + exp(1i*2*pi*f(k)*t)*reshape(conj(D(:,:,k)),1,[]);
    Dt = Dt + exp(-1i*2*pi*f(k)*t).*D(k,:,:) + exp(1i*2*pi*f(k)*t).*conj(D(k,:,:));
end

%%
imagesc(Dt(:,:,1) + Dt(:,:,100),[-1 1]*1e2)