% =========================================================================
% Aperiodic Fourier Modal Method
% =========================================================================
% Author: Ziwei Zhu
% =========================================================================

close all;
clear variables;

eta0 = -1i * 376.64;

% define materials
mu0 = 1.0;
eps0 = 2.25;
mu_r = 1.0;
eps_r = 12.25;

% simulation region
Cx = 0.0; 
Cy = 0.0; 
a = 1.00;
Lx = a; 
Ly = a;

% thickness
d = zeros(3, 1);
d(1) = 0.50;
d(2) = 0.50;
d(3) = 0.50;

% define source
freq = 3.001; % 2 * pi * c / a
lambda0 = a / (2 * pi * freq); % cm
k0  = 2.0 * pi / lambda0; % cm^-1

% number of harmonics
% M = ceil(7 * Lx / lambda0);
M = 11;
if mod(M, 2) == 0
    M = M + 1; 
end
    
%N = ceil(7 * Ly / lambda0);
N = 11;
if mod(N, 2) == 0
    N = N + 1;
end

ERC = zeros(M*N, M*N, 3);
URC = zeros(M*N, M*N, 3);

% grid parameters
Ny = 2 * max(M, N) + 1;
Nx = ceil((Lx / Ly) * Ny);
dx = Lx / Nx;
dy = Ly / Ny;

% build the geometry of the device
ER = ones(Nx, Ny, 3);
UR = ones(Nx, Ny, 3);
% %build device for layer 1
ER(:, :, 1) = eps0;
UR(:, :, 1) = mu0;
% %build device for layer 2
ER(:, :, 2) = eps0;
UR(:, :, 2) = mu0;
% %build device for layer 3
ER(:, :, 3) = eps0;
UR(:, :, 3) = mu0;

% [ER(:, :, 1), UR(:, :, 1)] =...
%     set_circle_pattern(Nx, Ny, Cx, Cy, Lx, Ly, 0.2, eps_r, eps0, mu_r, mu0);
% [ER(:, :, 2), UR(:, :, 2)] =...
%     set_circle_pattern(Nx, Ny, Cx, Cy, Lx, Ly, 0.2, eps_r, eps0, mu_r, mu0);
% [ER(:, :, 3), UR(:, :, 3)] =...
%     set_circle_pattern(Nx, Ny, Cx, Cy, Lx, Ly, 0.2, eps_r, eps0, mu_r, mu0);

ERC(:, :, 1) = convmat(ER(:, :, 1), M, N);
ERC(:, :, 2) = convmat(ER(:, :, 2), M, N);
ERC(:, :, 3) = convmat(ER(:, :, 3), M, N);
URC(:, :, 1) = convmat(UR(:, :, 1), M, N);
URC(:, :, 2) = convmat(UR(:, :, 2), M, N);
URC(:, :, 3) = convmat(UR(:, :, 3), M, N);

% build the geometry for waveguide excitation and acception
[ER_ref, UR_ref] = ...
    set_circle_pattern(Nx, Ny, Cx, Cy, Lx, Ly, 0.2, eps_r, eps0, mu_r, mu0);
[ER_trn, UR_trn] = ...
    set_circle_pattern(Nx, Ny, Cx, Cy, Lx, Ly, 0.2, eps_r, eps0, mu_r, mu0);

ERC_ref(:, :, 1) = convmat(ER_ref(:, :, 1), M, N);
ERC_trn(:, :, 1) = convmat(ER_trn(:, :, 1), M, N);

URC_ref(:, :, 1) = convmat(UR_ref(:, :, 1), M, N);
URC_trn(:, :, 1) = convmat(UR_trn(:, :, 1), M, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main RCWA routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup wave number matrices
wavenumber = struct;
wavenumber.M = floor(M / 2);
wavenumber.N = floor(N / 2);
wavenumber.m = -wavenumber.M : wavenumber.M;
wavenumber.n = -wavenumber.N : wavenumber.N;

% the code only works for normal incident
kxm = 0.0 - 2 * pi * wavenumber.m / k0 / Lx;
kyn = 0.0 - 2 * pi * wavenumber.n / k0 / Ly;

[Kyn, Kxm] = meshgrid(kyn, kxm);
KX = diag(Kxm(:));
KY = diag(Kyn(:));

% compute eigen modes of free space (analyze gap medium)
% because the thickness of the free space layer is zero, there should be no
% difference between the KX and KY with and without PML
% calculate W0, V0

W0 = eye(2*M*N);
Q0 = [KX * KY, eye(M*N) - KX * KX;
    KY * KY - eye(M*N), - KX * KY];
Kz0 = conj(sqrt(1.0 - Kxm.^2 - Kyn.^2));
KZ0 = diag(Kz0(:));
LAM0 = [1j * KZ0, zeros(M*N);
    zeros(M*N), 1j * KZ0];
V0 = Q0 / LAM0;

% multiply PML metric before KX and KY
c_PML = 1.0;
Dx = 0.9 * Lx;
Dy = 0.9 * Ly;

px = pml_coefficients(c_PML, Lx, Dx, 2*M + 1);
py = pml_coefficients(c_PML, Ly, Dy, 2*N + 1);

% PX = convmat_px(px, M, N, Lx);
% PY = convmat_py(py, M, N, Ly);
PX = eye(M*N);
PY = eye(M*N);

% initialize global s matrix
SG = struct;
SG.S11 = zeros(2*M*N);
SG.S12 = eye(2*M*N);
SG.S21 = eye(2*M*N);
SG.S22 = zeros(2*M*N);

% main loop through all layers
for i = 1 : 3
    % build eigen value problem
    ERi = ERC(:, :, i);
    URi = URC(:, :, i);
    
    Pi = [PX * KX / ERi * PY * KY, URi - PX * KX / ERi * PX * KX;
        PY * KY / ERi * PY * KY - URi, -PY * KY / ERi * PX * KX];
    Qi = [PX * KX / URi * PY * KY, ERi - PX * KX / URi * PX * KX;
        PY * KY / URi * PY * KY - ERi, -PY * KY / URi * PX * KX];
    
    OMEGA2 = Pi * Qi;
    
    [Wi, LAMi] = eig(OMEGA2);
    LAMi = sqrt_rectify(LAMi);
    Vi = Qi * Wi / LAMi;
    
    % calculate layer scatter matrix
    Ai0 = Wi \ W0 + Vi \ V0;
    Bi0 = Wi \ W0 - Vi \ V0;
    
    % should use expm instead of exp because exp(A) returns 1 for 0 entries
    % in A
    Xi = expm(-LAMi * k0 * d(i));
   
    S = struct;
    S.S11 = (Ai0 - Xi * Bi0 / Ai0 * Xi * Bi0)...
        \ (Xi * Bi0 / Ai0 * Xi * Ai0 - Bi0);
    S.S12 = (Ai0 - Xi * Bi0 / Ai0 * Xi * Bi0)...
        \ Xi * (Ai0 - Bi0 / Ai0 * Bi0);
    S.S21 = S.S12;
    S.S22 = S.S11;
    
    % update global scatter matrix
    SG = redheffer_star_product(SG, S);
    
end

% compute reflection side s matrix
Pref = [PX * KX / ERC_ref * PY * KY, URC_ref - PX * KX / ERC_ref * PX * KX;
    PY * KY / ERC_ref * PY * KY - URC_ref, -PY * KY / ERC_ref * PX * KX];
Qref = [PX * KX / URC_ref * PY * KY, ERC_ref - PX * KX / URC_ref * PX * KX;
    PY * KY / URC_ref * PY * KY - ERC_ref, -PY * KY / URC_ref * PX * KX];
OMEGAref = Pref * Qref;
[Wref, LAMref] = eig(OMEGAref);
LAMref = sqrt_rectify(LAMref);
Vref = Qref * Wref / LAMref;   

Aref = W0 \ Wref + V0 \ Vref;
Bref = W0 \ Wref - V0 \ Vref;

SR = struct;
SR.S11 = - Aref \ Bref;
SR.S12 = 2.0 * inv(Aref);
SR.S21 = 0.5 * (Aref - Bref / Aref * Bref);
SR.S22 = Bref / Aref;
SG = redheffer_star_product(SR, SG);

% compute transmission side s matrix
Ptrn = [PX * KX / ERC_trn * PY * KY, URC_trn - PX * KX / ERC_trn * PX * KX;
    PY * KY / ERC_trn * PY * KY - URC_trn, -PY * KY / ERC_trn * PX * KX];
Qtrn = [PX * KX / URC_trn * PY * KY, ERC_trn - PX * KX / URC_trn * PX * KX;
    PY * KY / URC_trn * PY * KY - ERC_trn, -PY * KY / URC_trn * PX * KX];
OMEGAtrn = Ptrn * Qtrn;
[Wtrn, LAMtrn] = eig(OMEGAtrn);
LAMtrn = sqrt_rectify(LAMtrn);
Vtrn = Qtrn * Wtrn / LAMtrn; 

Atrn = W0 \ Wtrn + V0 \ Vtrn;
Btrn = W0 \ Wtrn - V0 \ Vtrn;

ST = struct;
ST.S11 = Btrn / Atrn;
ST.S12 = 0.5 * (Atrn - Btrn / Atrn * Btrn);
ST.S21 = 2.0 * inv(Atrn);
ST.S22 = -Atrn \ Btrn;

% update global s matrix
SG = redheffer_star_product(SG, ST);

% solve for reflected & transmitted field

% excite from one polariz and for the fundamental mode
c_src = zeros(2*M*N, 1);
[~, center_mode_index] = max(abs(imag(diag(LAMref))));
% center_mode_index = 73;

% TEST
% for center_mode_index = 1 : 2*M*N
%     c_src = zeros(2*M*N, 1);
%     c_src(center_mode_index) = 1.0;
%     eT_src = Wref * c_src;
%     sx = eT_src(1:M*N);
%     sy = eT_src((M*N+1):2*M*N);
%     fig_temp1 = figure('Color', 'w');
%     imagesc(reshape(abs(sx), M, N));
%     fig_temp2 = figure('Color', 'w');
%     imagesc(reshape(abs(sy), M, N));
%     name1 = sprintf('modes_PML/sx_mode_%d.png', center_mode_index);
%     name2 = sprintf('modes_PML/sy_mode_%d.png', center_mode_index);
%     saveas(fig_temp1, name1);
%     saveas(fig_temp2, name2);
%     close all;
% end

c_src(center_mode_index) = 1.0;
c_ref = SG.S11 * c_src;
c_trn = SG.S21 * c_src;

eT_src = Wref * c_src;
ex_src = eT_src(1:M*N);
ey_src = eT_src((M*N+1):2*M*N);
hT_src = -eta0 * Vref * c_src;
hx_src = hT_src(1:M*N);
hy_src = hT_src((M*N+1):2*M*N);

eT_ref = Wref * c_ref;
ex_ref = eT_ref(1:M*N);
ey_ref = eT_ref((M*N+1):2*M*N);
hT_ref = eta0 * Vref * c_ref;
hx_ref = hT_ref(1:M*N);
hy_ref = hT_ref((M*N+1):2*M*N);

REF = -0.25 * ((sum(ex_ref .* conj(hy_ref) - ey_ref .* conj(hx_ref)))...
    + (sum(conj(ex_ref) .* hy_ref - conj(ey_ref) .* hx_ref))...
    + (sum(conj(ex_src) .* hy_ref - conj(ey_src) .* hx_ref))...
    + (sum(ex_ref .* conj(hy_src) - ey_ref .* conj(hx_src))));

SRC = Poynting_flux(c_src, Wref, Vref, eta0);
TRN = Poynting_flux(c_trn, Wtrn, Vtrn, eta0);

SUM = (-REF / SRC) + TRN / SRC;

% verify conservation
fprintf('SRC: %.3e, REF: %.3e, TRN: %.3e, SUM: %.3e\n', SRC, REF, TRN, SUM);


% here we only consider the fundamental mode
% c_ref(1 : (center_mode_index - 1)) = 0.0;
% c_ref((center_mode_index + 2) : 2*M*N) = 0.0;
% c_trn(1 : (center_mode_index - 1)) = 0.0;
% c_trn((center_mode_index + 2) : 2*M*N) = 0.0;

% % each column of Wref represents one mode, with one propagating constant
% eT_src = Wref * c_src;
% eT_ref = Wref * c_ref;
% eT_trn = Wtrn * c_trn;
% 
% % TODO: evaluate KZ_inc, KZ_ref, KZ_trn
% % we assume the wave are confined in the waveguide
% % this is true when we have a small amount of harmonics
% Kz_ref = -conj(sqrt(mu_r * eps_r - Kxm.^2 - Kyn.^2));
% Kz_trn = -Kz_ref;
% Kz_inc = Kz_trn;
% 
% KZ_ref = diag(Kz_ref(:));
% KZ_trn = diag(Kz_trn(:));
% KZ_inc = diag(Kz_inc(:));
% 
% sx = eT_src(1:M*N);
% sy = eT_src((M*N+1):2*M*N);
% sz = - KZ_inc \ (KX * sx + KY * sy);
% 
% rx = eT_ref(1:M*N);
% ry = eT_ref((M*N+1):2*M*N);
% rz = - KZ_ref \ (KX * rx + KY * ry);
% 
% tx = eT_trn(1:M*N);
% ty = eT_trn((M*N+1):2*M*N);
% tz = - KZ_trn \ (KX * tx + KY * ty);
% 
% S2 = abs(sx).^2 + abs(sy).^2 + abs(sz).^2;
% Senergy = real(KZ_inc / mu0) * S2;
% Senergy = reshape(Senergy, M, N);
% SRC = sum(Senergy(:));
% 
% % calculate diffraction efficiencies
% R2 = abs(rx).^2 + abs(ry).^2 + abs(rz).^2;
% R = real(-KZ_ref / mu0)  * R2;
% R = reshape(R, M, N);
% REF = sum(R(:)) / SRC;
% 
% T2 = abs(tx).^2 + abs(ty).^2 + abs(tz).^2;
% T = real(KZ_trn / mu0) * T2;
% T = reshape(T, M, N);
% TRN = sum(T(:)) / SRC;
% 
% SUM = REF + TRN;
% 
% % verify conservation
% fprintf('REF: %.3f, TRN: %.3f, SUM: %.3f\n', REF, TRN, SUM);

%end

% fig1 = figure('Color', 'w');
% imagesc(Senergy);
% colorbar;
% axis equal tight;
% saveas(fig1, 'source.png');
% 
% fig2 = figure('Color', 'w');
% imagesc(R);
% colorbar;
% axis equal tight;
% saveas(fig2, 'reflectance.png');
% 
% fig3 = figure('Color', 'w');
% imagesc(T);
% colorbar;
% axis equal tight;
% saveas(fig3, 'Transmittance.png');
