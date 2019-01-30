%% setup

% computational domain
dx = 10;
L = [2000, 3000];

% define model
v = @(z,x) 2 + 0.5*(z>500) + 0.5*(z + 100*cos(1e-2*x) > 800) + 0.5*(z - 0.1*x > 1500);

% initial model
v0 = @(z,x)2 + 1e-3*z;

% set frequency, wavelet etc.
t = 0:1e-3:2;
nt = length(t);
f = 0:1/t(end):.5/(t(2)-t(1));
f0 = 10;
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
nsrc = 50;
zs = .1*L(1)*ones(1,nsrc);
xs = linspace(.25*L(2),.75*L(2),nsrc);

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

%D = F(m,Q,model);
%save('data.mat',D);

load('data.mat')

%% Blending
nsim = 10;

phi = 2*pi*f(If)*(.1 + 0*rand(length(If)));
%Bs = Q(mod(0 + (0:5:nsrc-1),nsrc)+1,:) + Q(mod(1 + (0:5:nsrc-1),nsrc)+1,:) + Q(mod(2 + (0:5:nsrc-1),nsrc)+1,:) + Q(mod(3 + (0:5:nsrc-1),nsrc)+1,:) + Q(mod(4 + (0:5:nsrc-1),nsrc)+1,:);
Bs = randn(nsim,nrec);
Br = opDirac(nrec);
Bf = opDiag(exp(1i*phi));

B = opKron(opMatrix(Bs),Br);

Y = B*D*Bf';

%% deblending
X = zeros(nrec*nsrc,length(If));
for k = 1:length(nf)
    Yk = reshape(Y(:,k),nrec,nsim);
    Xk = exp(1i*phi(k))*Br'*Yk*Bs;
    X(:,k) = Xk(:);
end

%X = B'*Y*Bf;

%% Fourier transform
Fr = opDFTR(nt);
Rf = opRestriction(length(f),If);
D = reshape(D,nsrc*nrec,[]);
W = opDiag(Fr*wt(:));

Dt = Fr'*W*Rf'*D';
Yt = Fr'*W*Rf'*Y';
Xt = Fr'*W*Rf'*X';

%%
Is = [1 nsrc/2 nsrc];
Iss = [1 nsim/2 nsim];
Dt = reshape(Dt,nt,nrec,nsrc);
Yt = reshape(Yt,nt,nrec,nsim);
Xt = reshape(Xt,nt,nrec,nsrc);

figure;
imagesc(x,z,v(zz,xx));axis equal tight

figure;
Dt = reshape(Dt,nt,nrec,nsrc);
for k = 1:3
    subplot(1,3,k);imagesc(xr,t,Dt(:,:,Is(k)));colormap(gray)
end

figure;
Yt = reshape(Yt,nt,nrec,nsim);
for k = 1:3
    subplot(1,3,k);imagesc(xr,t,Yt(:,:,Iss(k)));colormap(gray)
end

figure;
Xt = reshape(Xt,nt,nrec,nsrc);
for k = 1:3
    subplot(1,3,k);imagesc(xr,t,Xt(:,:,Is(k)));colormap(gray)
end
%%