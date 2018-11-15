% =========================================================================
% 2D Finite Difference Beam Propagation Method (E Mode)
% =========================================================================
% Author: Ziwei Zhu
% Reference: http://emlab.utep.edu/ee5390cem.htm
% =========================================================================
close all;
clear variables;

% physical constant

% waveguide inner index (SiO2)
n1 = 1.46;
% waveguide outer index (Si)
n2 = 3.48;

% effective propagation refractive index
n_eff = 2.0;

% wavelength of light in vacuum
lambda0 = 1.55e-6; 

k0 = 2 * pi / lambda0;

% define simulation parameters
z0 = 0.0;
z1 = 2e-4;
x0 = -1e-5;
x1 = 1e-5;

wg_1_z0 = 0.0;
wg_1_z1 = 1.5e-4;
wg_1_x0 = -5e-6;
wg_1_x1 = -1e-6;

wg_2_z0 = 5e-5;
wg_2_z1 = 2e-4;
wg_2_x0 = 1e-6;
wg_2_x1 = 5e-6;


% normalize w.r.t wave vector
z0 = z0 * k0;
z1 = z1 * k0;
x0 = x0 * k0;
x1 = x1 * k0;

wg_1_z0 = wg_1_z0 * k0;
wg_1_z1 = wg_1_z1 * k0;
wg_1_x0 = wg_1_x0 * k0;
wg_1_x1 = wg_1_x1 * k0;

wg_2_z0 = wg_2_z0 * k0;
wg_2_z1 = wg_2_z1 * k0;
wg_2_x0 = wg_2_x0 * k0;
wg_2_x1 = wg_2_x1 * k0;

grid_per_devlen = 16;
grid_per_wavelen = 20;

dev_min_feature_len = min(abs(wg_1_x1 - wg_1_x0), abs(wg_2_x1 - wg_2_x0));

res_wave = 2 * pi / grid_per_wavelen / n_eff;
res_dev = dev_min_feature_len / grid_per_devlen;

Nx = ceil((x1 - x0) / res_dev);
Nz = ceil((z1 - z0) / res_wave);

res_z = (z1 - z0) / Nz;
res_x = (x1 - x0) / Nx;


% build device on grid
Mu_xx = ones(Nz, Nx);
Mu_zz = ones(Nz, Nx);

Eps_yy = ones(Nz, Nx);

eps1 = n1 * n1;
eps2 = n2 * n2;

for i = 1 : Nz
    for j = 1 : Nx
        z = z0 + res_z * (i - 0.5);
        x = x0 + res_x * (j - 0.5);
        
        if (z >= wg_1_z0 && z <= wg_1_z1 && x <= wg_1_x1 && x >= wg_1_x0)
            Eps_yy(i, j) = eps1;
        elseif (z >= wg_2_z0 && z <= wg_2_z1 && x <= wg_2_x1 && x >= wg_2_x0)
            Eps_yy(i, j) = eps1;
        else
            Eps_yy(i, j) = eps2;
        end
    end
end

% TODO handle PML layers


% compute matrix derivative operators
% use central difference
% assume fully reflectance boundary
Dx_h = (0.5 * diag(ones(Nx-1,1), 1) - 0.5 * diag(ones(Nx-1, 1), -1)) / res_x;
Dx_e = (0.5 * diag(ones(Nx-1,1), 1) - 0.5 * diag(ones(Nx-1, 1), -1)) / res_x;


% initialize field
Ey = zeros(Nz, Nx);
Hx = zeros(Nz, Nx);
Hz = zeros(Nz, Nx);

% compute source at the first plane
tau = 0.1 * (wg_1_x1 - wg_1_x0);
delay_pos =  0.5 * (wg_1_x1 + wg_1_x0);
Generator = @(t)exp(-0.5 * ((t - delay_pos) ./ tau).^2);
for i = 1 : Nx
    x = x0 + res_x * (i - 0.5);
    Ey(1, i) = Generator(x);
end

% iterate to update field
for T = 1 : Nz
    
end