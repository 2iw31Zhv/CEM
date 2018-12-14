% =========================================================================
% Rigorous Coupled Wave Analysis
% =========================================================================
% Author: Ziwei Zhu
% Reference: http://emlab.utep.edu/ee5390cem.htm
% =========================================================================

% define materials
mu_ref = 1.0;
eps_ref = 2.0;
mu_r = 1.0;
eps_r = 6.0;
mu_trn = 1.0;
eps_trn = 9.0;

% define geometries

% simulation region
Cx = 0.0; % cm
Cy = 0.0; % cm
Lx = 1.75; % cm
Ly = 1.50; % cm

% thickness
d = zeros(2, 1);
d(1) = 0.50; % cm
d(2) = 0.30; % cm

% hole size of layer 1
w = 0.8 * Lx;


% define source
lambda0 = 2.0; % cm
k0  = 2.0 * pi / lambda0; % cm^-1
theta_src = 0; % arc
phi_src = pi / 4.0; % arc


n_inc = sqrt(eps_ref / mu_ref);
k_inc = n_inc * [ sin(theta_src) * cos(phi_src);
    sin(theta_src) * sin(phi_src);
    cos(theta_src)];


% grid parameters
Ny = 16;
Nx = ceil((Lx / Ly) * Ny);
dx = Lx / Nx;
dy = Ly / Ny;

% build the geometry of the device
ER = ones(Nx, Ny, 2);
UR = ones(Nx, Ny, 2);

% build device for layer 1
for i  = 1 : Nx
    for j = 1 : Ny
        x = Cx - 0.5 * Lx + (i - 0.5) * dx;
        y = Cy - 0.5 * Ly + (j - 0.5) * dy;
        
        if (y <= sqrt(3) * x + sqrt(3) / 4 * w && ...
            y <= -sqrt(3) * x + sqrt(3) / 4 * w && ...
            y >= -sqrt(3) / 4 * w)
            ER(i, j, 1) = eps_ref;
            UR(i, j, 1) = mu_ref;
        else
            ER(i, j, 1) = eps_r;
            UR(i, j, 1) = mu_r;
        end
    end
end

% build device for layer 2
ER(:, :, 2) = eps_r;
UR(:, :, 2) = mu_r;

% calculate the corresponding convolution matrix;

% number of harmonics
M = ceil(7 * Lx / lambda0);
if mod(M, 2) == 0
    M = M + 1; 
end
    
N = ceil(7 * Ly / lambda0);
if mod(N, 2) == 0
    N = N + 1;
end

ERC = zeros(M*N, M*N, 2);
URC = zeros(M*N, M*N, 2);

ERC(:, :, 1) = convmat(ER(:, :, 1), M, N);
ERC(:, :, 2) = convmat(ER(:, :, 2), M, N);

URC(:, :, 1) = convmat(UR(:, :, 1), M, N);
URC(:, :, 2) = convmat(UR(:, :, 2), M, N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main RCWA routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup wave number matrices
wavenumber = struct;
wavenumber.M = floor(M / 2);
wavenumber.N = floor(N / 2);
wavenumber.m = -wavenumber.M : wavenumber.M;
wavenumber.n = -wavenumber.N : wavenumber.N;

kxm = k_inc(1) - 2 * pi * wavenumber.m / k0 / Lx;
kyn = k_inc(2) - 2 * pi * wavenumber.n / k0 / Ly;

[Kyn, Kxm] = meshgrid(kyn, kxm);

Kz_ref = -sqrt(mu_ref' * eps_ref' - Kxm.^2 - Kyn.^2)';
Kz_trn = sqrt(mu_trn' * eps_trn' - Kxm.^2 - Kyn.^2)';

KX = diag(sparse(Kxm(:)));
KY = diag(sparse(Kyn(:)));
KZ_ref = diag(sparse(Kz_ref(:)));
KZ_trn = diag(sparse(Kz_trn(:)));

% compute eigen modes of free space (analyze gap medium)
% calculate W0, V0

W0 = sparse(eye(2*M*N));
Q0 = [KX * KY, eye(M*N) - KX * KX;
    KY * KY - eye(M*N), - KX * KY];
Kz0 = sqrt(1.0 - Kxm.^2 - Kyn.^2)';
KZ0 = diag(sparse(Kz0(:)));
LAM0 = [1j * KZ0, zeros(M*N);
    zeros(M*N), 1j * KZ0];
V0 = Q0 / LAM0;


% initialize global s matrix
SG = struct;
SG.S11 = zeros(2*M*N);
SG.S12 = eye(2*M*N);
SG.S21 = eye(2*M*N);
SG.S22 = zeros(2*M*N);

% main loop through all layers
for i = 1 : 2
    % build eigen value problem
    ERi = ERC(:, :, i);
    URi = URC(:, :, i);
    
    Pi = [KX / ERi * KY, URi - KX / ERi * KX;
        KY / ERi * KY - URi, -KY / ERi * KX];
    Qi = [KX / URi * KY, ERi - KX / URi * KX;
        KY / URi * KY - ERi, -KY / URi * KX];
    OMEGA2 = Pi * Qi;
    
    [Wi, LAMi] = eig(OMEGA2);
    LAMi = sqrt(LAMi);
    Vi = Qi * Wi / LAMi;
    
    % calculate layer scatter matrix
    Ai0 = Wi \ W0 + Vi \ V0;
    Bi0 = Wi \ W0 - Vi \ V0;
    Xi = exp(-LAMi * k0 * d(i));
    
    
    S = struct;
    S.S11 = (Ai0 - Xi * Bi0 / Ai0 * Xi * Bi0)...
        \ (Xi * (Bi0 / Ai0) * Xi * Ai0 - Bi0);
    S.S12 = (Ai0 - Xi * Bi0 / Ai0 * Xi * Bi0)...
        \ Xi * (Ai0 - Bi0 / Ai0 * Bi0);
    S.S21 = S.S12;
    S.S22 = S.S11;
    
    % update global scatter matrix
    SG = redheffer_star_product_expand2(SG, S);
    
end

% compute reflection side s matrix
Qref = 1.0 / mu_ref * [KX * KY, mu_ref * eps_ref * eye(M*N) - KX * KX;
    KY * KY - mu_ref * eps_ref * eye(M*N), -KY * KX];
Wref = eye(2*M*N);
LAMref = [-1j * KZ_ref, zeros(M*N);
    zeros(M*N), -1j * KZ_ref];
Vref = Qref / LAMref;

Aref = W0 \ Wref + V0 \ Vref;
Bref = W0 \ Wref - V0 \ Vref;

SR = struct;
SR.S11 = - Aref \ Bref;
SR.S12 = 2.0 * inv(Aref);
SR.S21 = 0.5 * (Aref - Bref / Aref * Bref);
SR.S22 = Bref / Aref;

% compute transmission side s matrix
Qtrn = 1.0 / mu_trn * [KX * KY, mu_trn * eps_trn * eye(M*N) - KX * KX;
    KY * KY - mu_trn * eps_trn * eye(M*N), -KY * KX];
Wtrn = eye(2*M*N);
LAMtrn = [1j * KZ_trn, zeros(M*N);
    zeros(M*N), 1j * KZ_trn];
Vtrn = Qtrn / LAMtrn;

Atrn = W0 \ Wtrn + V0 \ Vtrn;
Btrn = W0 \ Wtrn - V0 \ Vtrn;

ST = struct;
ST.S11 = Btrn / Atrn;
ST.S12 = 0.5 * (Atrn - Btrn / Atrn * Btrn);
ST.S21 = 2.0 * inv(Atrn);
ST.S22 = -Atrn \ Btrn;

% update global s matrix
SG = redheffer_star_product_expand2(SR, SG);
SG = redheffer_star_product_expand2(SG, ST);

% compute source
source = struct;

source.delta = zeros(M*N, 1);
source.delta(floor(N/2) * M + floor(M/2) + 1) = 1.0;

source.n = [0; 0; 1];
% for normal incidence
source.vTE = [0; 1; 0];
source.vTM = cross(k_inc, source.vTE);
source.vTM = source.vTM / norm(source.vTM);

source.theta = pi / 6.0;
source.P = cos(source.theta) * source.vTM + sin(source.theta) * source.vTE;

source.eT = [source.P(1) * source.delta;
    source.P(2) * source.delta];

% solve for reflected & transmitted field
c_src = Wref \ source.eT;
c_ref = SG.S11 * c_src;
c_trn = SG.S21 * c_src;

eT_ref = Wref * c_ref;
eT_trn = Wtrn * c_trn;

rx = eT_ref(1:M*N);
ry = eT_ref((M*N+1):2*M*N);
rz = - KZ_ref \ (KX * rx + KY * ry);

tx = eT_trn(1:M*N);
ty = eT_trn((M*N+1):2*M*N);
tz = - KZ_trn \ (KX * tx + KY * ty);

% calculate diffraction efficiencies
R2 = abs(rx).^2 + abs(ry).^2 + abs(rz).^2;
R = real(-KZ_ref / mu_ref) / real(k_inc(3) / mu_ref) * R2;
R = reshape(R, M, N);

REF = sum(R(:));


T2 = abs(tx).^2 + abs(ty).^2 + abs(tz).^2;
T = real(KZ_trn / mu_trn) / real(k_inc(3) / mu_ref) * T2;
T = reshape(T, M, N);

TRN = sum(T(:));

% verify conservation

REF + TRN
