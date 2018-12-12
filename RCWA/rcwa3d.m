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


% physics of the field
lambda0 = 2.0; % cm

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
Q = ceil(7 * Ly / lambda0);

ERC = zeros(P*Q, P*Q, 2);
URC = zeros(P*Q, P*Q, 2);

ERC(:, :, 1) = convmat(ER(:, :, 1), P, Q);
ERC(:, :, 2) = convmat(ER(:, :, 2), P, Q);

URC(:, :, 1) = convmat(UR(:, :, 1), P, Q);
URC(:, :, 2) = convmat(UR(:, :, 2), P, Q);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main RCWA routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







