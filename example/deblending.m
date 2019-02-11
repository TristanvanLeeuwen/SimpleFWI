%% setup

% computational domain
dx = 50;
L = [2000, 3000];

% define model
v = @(z,x) 2 + 0.5*(z>500) + 0.5*(z + 100*cos(1e-2*x) > 800) + 0.5*(z - 0.1*x > 1500);

% initial model
v0 = @(z,x)2 + 1e-3*z;

% set frequency, wavelet etc.
t = 0:1e-3:2;
nt = length(t);
f = 0:1/t(end):.5/(t(2)-t(1));
f0 = 5;
t0 = 0.1;
wt = (1 - 2*(pi*f0*(t-t0)).^2).*exp(-(pi*f0*(t-t0)).^2);
wf = (2/sqrt(pi)/f0).*(f/f0).^2.*exp(-(f/f0).^2);
If = find(wf/max(wf)>1e-3);
nf = length(If);

% receivers
nrec = 50;
zr = .1*L(1)*ones(1,nrec);
xr = linspace(.25*L(2),.75*L(2),nrec);


% sources
nsrc = 2;
zs = .1*L(1)*ones(1,nsrc);
xs = linspace(.4*L(2),.6*L(2),nsrc);

% regularization parameter
alpha = 0;

%% observed data
% grid
n  = L./dx + 1;
h  = dx*[1 1];
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

% parameters
model.f = f(If);
model.w = wf(If);
model.n = n;
model.h = h;
model.zr = zr;
model.xr = xr;
model.zs = zs;
model.xs = xs;

% model
m = 1./v(zz(:),xx(:)).^2;

% data
Q = eye(nsrc);

D = F(m,Q,model);
save('data.mat','D');

%load('data.mat')

%% Blending
nsim = 2;
a = 1 + 2*rand(nsrc,nsim);
tau = .1 + .2*rand(nsrc,nsim);

B = @(f) Q*a.*exp(1i*2*pi*f*tau);

D = reshape(D,nrec,nsrc,length(If));
Y = zeros(nrec,nsim,length(If));
for k = 1:length(If)
    Y(:,:,k) = D(:,:,k)*B(f(If(k)));
end


%% deblending
Y = reshape(Y,nrec,nsim,length(If));
X = zeros(nrec, nsrc,length(If));
for k = 1:length(If)
    %X(:,:,k) = Y(:,:,k)*B(f(If(k)))';
    X(:,:,k) = Y(:,:,k)*B(f(If(k)))'*inv(B(f(If(k)))*B(f(If(k)))');
end

%% Fourier transform
Fr = opDFTR(nt);
Rf = opRestriction(length(f),If);

D = reshape(D,nsrc*nrec,[]);
Y = reshape(Y,nsim*nrec,[]);
X = reshape(X,nsrc*nrec,[]);

Dt = Fr'*Rf'*D';
Yt = Fr'*Rf'*Y';
Xt = Fr'*Rf'*X';

%%
Is = [1,2];
Iss = [1,2];
Dt = reshape(Dt,nt,nrec,nsrc);
Yt = reshape(Yt,nt,nrec,nsim);
Xt = reshape(Xt,nt,nrec,nsrc);

figure;
imagesc(x,z,v(zz,xx));axis equal tight

figure;
Dt = reshape(Dt,nt,nrec,nsrc);
for k = 1:length(Is)
    subplot(1,length(Is),k);imagesc(xr,t,Dt(:,:,Is(k)));colormap(gray)
end

figure;
Yt = reshape(Yt,nt,nrec,nsim);
for k = 1:length(Iss)
    subplot(1,length(Iss),k);imagesc(xr,t,Yt(:,:,Iss(k)));colormap(gray)
end

figure;
Xt = reshape(Xt,nt,nrec,nsrc);
for k = 1:length(Is)
    subplot(1,length(Is),k);imagesc(xr,t,Xt(:,:,Is(k)));colormap(gray)
end
%%