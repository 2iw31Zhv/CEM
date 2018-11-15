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
n_eff = 2.8;

% wavelength of light in vacuum
lambda0 = 1.55e-6; 

k0 = 2 * pi / lambda0;

% define simulation parameters
z0 = 0.0;
z1 = 2e-5;
x0 = -1e-5;
x1 = 1e-5;

wg_1_z0 = 0.0;
wg_1_z1 = 1.5e-5;
wg_1_x0 = -5e-6;
wg_1_x1 = -1e-6;

wg_2_z0 = 5e-6;
wg_2_z1 = 2e-5;
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

grid_per_devlen = 32;
grid_per_wavelen = 100;

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

% compute source at the first plane
tau = 0.5 * (wg_1_x1 - wg_1_x0);
delay_pos =  0.5 * (wg_1_x1 + wg_1_x0);
Generator = @(t)exp(-0.5 * ((t - delay_pos) ./ tau).^2);
for i = 1 : Nx
    x = x0 + res_x * (i - 0.5);
    Ey(1, i) = Generator(x);
end

% iterate to update field
for i = 1 : Nz-1
    % extract and diagonalize constitutive tensor
    mu_xx0 = diag(Mu_xx(i, :));
    mu_zz0 = diag(Mu_zz(i, :));
    eps_yy0 = diag(Eps_yy(i, :));
    
    inverse_mu_zz0 = diag(Mu_zz(i,:).^(-1));
    
    A0 = mu_xx0 * Dx_h * inverse_mu_zz0 * Dx_e + mu_xx0 * eps_yy0...
        - n_eff^2 * ones(Nx, Nx);
    
    mu_xx1 = diag(Mu_xx(i+1, :));
    mu_zz1 = diag(Mu_zz(i+1, :));
    eps_yy1 = diag(Eps_yy(i+1, :));
    
    inverse_mu_zz1 = diag(Mu_zz(i+1,:).^(-1));
    
    A1 = mu_xx1 * Dx_h * inverse_mu_zz1 * Dx_e + mu_xx1 * eps_yy1...
        - n_eff^2 * ones(Nx, Nx);  
    
    left_hand_mat = ones(Nx, Nx) - 1j * res_z / 4.0 / n_eff * A1;
    rhs = (ones(Nx, Nx) + 1j * res_z / 4.0 / n_eff * A0) * (Ey(i, :)');
    
    Ey(i+1, :) = (left_hand_mat \ rhs)';
    x_array = x0 + res_x * ((1:Nx)-0.5);
    z_array = z0 + res_z * ((1:Nz)-0.5);
    
    if (mod(i, 20) == 0)
        clear gca;
        imagesc(z_array, x_array, abs(Ey)');
        hold on;
        plot([wg_1_z0 wg_1_z0], [wg_1_x0 wg_1_x1], '-k');
        hold on;
        plot([wg_1_z0 wg_1_z1], [wg_1_x0 wg_1_x0], '-k');
        hold on;
        plot([wg_1_z1 wg_1_z1], [wg_1_x0 wg_1_x1], '-k');
        hold on;
        plot([wg_1_z0 wg_1_z1], [wg_1_x1 wg_1_x1], '-k');
        hold on;
        plot([wg_2_z0 wg_2_z0], [wg_2_x0 wg_2_x1], '-k');
        hold on;
        plot([wg_2_z0 wg_2_z1], [wg_2_x0 wg_2_x0], '-k');
        hold on;
        plot([wg_2_z1 wg_2_z1], [wg_2_x0 wg_2_x1], '-k');
        hold on;
        plot([wg_2_z0 wg_2_z1], [wg_2_x1 wg_2_x1], '-k');
        pause(0.01);
    end
end