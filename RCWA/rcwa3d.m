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
d1 = 0.50; % cm
d2 = 0.30; % cm

% hole size of layer 1
w = 0.8 * Lx;


% define source
lambda0 = 2.0; % cm
k0  = 2.0 * pi / lambda0; % cm^-1
theta_src = pi / 4.0; % arc
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

imagesc(ER(:, :, 1)');
axis equal tight;

% build device for layer 2
ER(:, :, 2) = eps_r;
UR(:, :, 2) = mu_r;

% calculate the corresponding convolution matrix;

% number of harmonics
P = ceil(7 * Lx / lambda0);
if mod(P, 2) == 0
    P = P + 1; 
end
    
Q = ceil(7 * Ly / lambda0);
if mod(Q, 2) == 0
    Q = Q + 1;
end

ERC = zeros(P*Q, P*Q, 2);
URC = zeros(P*Q, P*Q, 2);

ERC(:, :, 1) = convmat(ER(:, :, 1), P, Q);
ERC(:, :, 2) = convmat(ER(:, :, 2), P, Q);

URC(:, :, 1) = convmat(UR(:, :, 1), P, Q);
URC(:, :, 2) = convmat(UR(:, :, 2), P, Q);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main RCWA routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup wave number matrices
wavenumber = struct;
wavenumber.M = floor(P / 2);
wavenumber.N = floor(Q / 2);
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

% compute eigen modes of free space




