% =========================================================================
% 2D Finite Difference Beam Propagation Method (E Mode)
% =========================================================================
% Author: Ziwei Zhu
% Reference: http://emlab.utep.edu/ee5390cem.htm
% =========================================================================
close all;
clear variables;

% physical constant

% waveguide inner index (Si)
n1 = 3.48;
% waveguide outer index (SiO2)
n2 = 1.46;

% wavelength of light in vacuum
lambda0 = 1.55e-6; 

k0 = 2 * pi / lambda0;

% define simulation parameters
z0 = 0.0;
z1 = 2e-5;
x0 = -1e-5;
x1 = 1e-5;

wg_1_z0 = 0.0;
wg_1_z1 = 2e-5;
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


% effective propagation refractive index
% for the first mode:
sin_theta = binary_search_sin_theta(n1, n2, abs(wg_1_x1 - wg_1_x0), 2*pi);
n_eff = n1 * sqrt(1 - sin_theta^2);

grid_per_devlen = 4;
grid_per_wavelen = 10;

dev_min_feature_len = min(abs(wg_1_x1 - wg_1_x0), abs(wg_2_x1 - wg_2_x0));

res_wave = 2 * pi / grid_per_wavelen / n_eff;
res_dev = dev_min_feature_len / grid_per_devlen;
res = min(res_wave, res_dev);

Nx = ceil((x1 - x0) / res);
Nz = ceil((z1 - z0) / (0.1 *res));

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

% handle PML layers
sx = ones(Nz, Nx);
NPML_x = 40;
a_max = 3.0;
PML_power = 3;
eta0 = 376.73;

for i = 1 : NPML_x
    nx = NPML_x - i + 1;
    sx(:, i) = (1 + a_max * (nx / NPML_x)^PML_power)...
        * (1 + 1j * eta0 * (sin(pi * nx / 2 / NPML_x))^2);
end

for i = (Nx - NPML_x+1) : Nx
    nx = NPML_x - (Nx - i);
    sx(:, i) = (1 + a_max * (nx / NPML_x)^PML_power)...
        * (1 + 1j * eta0 * (sin(pi * nx / 2 / NPML_x))^2);
end
Eps_yy = Eps_yy .* sx;
Mu_xx = Mu_xx ./ sx;
Mu_zz = Mu_zz .* sx;

% compute matrix derivative operators
% use central difference
% assume fully reflectance boundary

% backward difference
Dx_e = (diag(ones(Nx,1)) - diag(ones(Nx-1, 1), -1)) / res_x;
% this is for the periodic boundary condition
% Dx_e(1, Nx) = -1 / res_x;

% forward difference
Dx_h = (-diag(ones(Nx,1)) + diag(ones(Nx-1, 1), 1)) / res_x;
% this is for the periodic boundary condition
% Dx_h(Nx, 1) = 1 / res_x;

% initialize field
Ey = zeros(Nz, Nx);

% compute source at the first plane
tau = 0.5 * (wg_1_x1 - wg_1_x0);
delay_pos =  0.5 * (wg_1_x1 + wg_1_x0);
Generator = @(t)exp(-0.5 * ((t - delay_pos) ./ tau).^2);

for i = 1 : Nx
    x = x0 + res_x * (i - 0.5);
    %if (x >= wg_1_x0 && x <= wg_1_x1)
        Ey(1, i) = Generator(x);
    %end
end

fig = figure('Color', 'w', 'Position', [300 300 500 400]);
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
    [Z, X] = meshgrid(z_array, x_array);
    phase = exp(1i * n_eff * Z);
    
    if (mod(i, 20) == 0)
        clear gca;
        
        h(1) = subplot(1, 2, 1);
        plot(x_array, Ey(1, :));
        camroll(-90);
        axis tight;
        
        h(2) = subplot(1, 2, 2);
        imagesc(z_array, x_array, abs(Ey)');
        axis equal tight;
        caxis([0, 1]);
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
        
        pos0 = get(h(1), 'Position');
        pos = get(h(2), 'Position');
        set(h(1), 'Position', [0.1, 0.15, 0.1, 0.75]);
        set(h(2), 'Position',[0.2, 0.15, 0.75, 0.75]);
        
        pause(0.01);
    end
end

saveas(gcf, 'beam_prop.png')