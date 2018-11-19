% =========================================================================
% Generator of Example Waveguide Coupler (E Mode)
% =========================================================================
% Author: Ziwei Zhu
% =========================================================================
close all;
clear variables;

% waveguide inner index (Si)
n1 = 3.48;
% waveguide outer index (SiO2)
n2 = 1.46;

% definition of geometry
width_1 = 1.95e-6;
zc_ellipse_1 = 5e-6;
xc_ellipse_1 = -1e-5;

a1_ellipse_1 = 1.2e-5;
a0_ellipse_1 = a1_ellipse_1 - width_1;

b1_ellipse_1 = 1e-5 + 0.5 * width_1;
b0_ellipse_1 = 1e-5 - 0.5 * width_1;

width_2 = 1.95e-6;
zc_ellipse_2 = 1.5e-5;
xc_ellipse_2 = 1e-5;

a1_ellipse_2 = 1.2e-5;
a0_ellipse_2 = a1_ellipse_2 - width_2;

b1_ellipse_2 = 1e-5 + 0.5 * width_2;
b0_ellipse_2 = 1e-5 - 0.5 * width_2;


% simulation defintion
Nz = 408;
Nx = 408;

z0 = 0.0;
z1 = 2e-5;
x0 = -1e-5;
x1 = 1e-5;

res_z = (z1 - z0) / Nz;
res_x = (x1 - x0) / Nx;


eps1 = n1 * n1;
eps2 = n2 * n2;

Mu_xx = ones(Nz, Nx);
Mu_zz = ones(Nz, Nx);
Eps_yy = ones(Nz, Nx);

for i = 1 : Nz
    for j = 1 : Nx
        z = z0 + res_z * (i - 0.5);
        x = x0 + res_x * (j - 0.5);
        if (((z - zc_ellipse_1)^2 / a1_ellipse_1^2 ...
            + (x - xc_ellipse_1)^2 / b1_ellipse_1^2 <= 1.0)...
            &&...
            ((z - zc_ellipse_1)^2 / a0_ellipse_1^2 ...
            + (x - xc_ellipse_1)^2 / b0_ellipse_1^2 >= 1.0)...
            &&...
            z >= zc_ellipse_1)
            Eps_yy(i, j) = eps1;
        elseif (z <= zc_ellipse_1 && x <= 0.5 * width_1 && x >= -0.5 * width_1)
            Eps_yy(i, j) = eps1;
        elseif (((z - zc_ellipse_2)^2 / a1_ellipse_2^2 ...
            + (x - xc_ellipse_2)^2 / b1_ellipse_2^2 <= 1.0)...
            &&...
            ((z - zc_ellipse_2)^2 / a0_ellipse_2^2 ...
            + (x - xc_ellipse_2)^2 / b0_ellipse_2^2 >= 1.0)...
            &&...
            z <= zc_ellipse_2)
            Eps_yy(i, j) = eps1;
        elseif (z >= zc_ellipse_2 && x <= 0.5 * width_1 && x >= -0.5 * width_1)
            Eps_yy(i, j) = eps1;
        else
            Eps_yy(i, j) = eps2;    
        end
    end
    
end

N_plot = 100;

wg_1_z0 = zeros(N_plot, 1);
wg_1_x0 = zeros(N_plot, 1);
wg_1_z1 = zeros(N_plot, 1);
wg_1_x1 = zeros(N_plot, 1);

wg_2_z0 = zeros(N_plot, 1);
wg_2_x0 = zeros(N_plot, 1);
wg_2_z1 = zeros(N_plot, 1);
wg_2_x1 = zeros(N_plot, 1);

for i = 1 : ceil(0.25*N_plot)-1
    alpha = i / (ceil(0.25*N_plot)-1);
    wg_1_z0(i) = alpha * zc_ellipse_1 + (1 - alpha) * z0;
    wg_1_x0(i) = -0.5 * width_1;
    wg_1_z1(i) = alpha * zc_ellipse_1 + (1 - alpha) * z0;
    wg_1_x1(i) = 0.5 * width_1;
    
    wg_2_z0(i) = alpha * zc_ellipse_2 + (1 - alpha) * z1;
    wg_2_x0(i) = 0.5 * width_1;
    wg_2_z1(i) = alpha * zc_ellipse_2 + (1 - alpha) * z1;
    wg_2_x1(i) = -0.5 * width_1;
end

for i = ceil(0.25 * N_plot) : N_plot
    alpha = (i - ceil(0.25 * N_plot) + 1) / (N_plot - ceil(0.25 * N_plot)+1);
    wg_1_z0(i) = alpha * (zc_ellipse_1 + a0_ellipse_1)...
        + (1 - alpha) * zc_ellipse_1;
    wg_1_x0(i) = xc_ellipse_1...
        + b0_ellipse_1 * sqrt(1 - (wg_1_z0(i) - zc_ellipse_1)^2 / a0_ellipse_1^2);
    wg_1_z1(i) = alpha * (zc_ellipse_1 + a1_ellipse_1)...
        + (1 - alpha) * zc_ellipse_1;
    wg_1_x1(i) = xc_ellipse_1...
        + b1_ellipse_1 * sqrt(1 - (wg_1_z1(i) - zc_ellipse_1)^2 / a1_ellipse_1^2);
    
    wg_2_z0(i) = alpha * (zc_ellipse_2 - a0_ellipse_2)...
        + (1 - alpha) * zc_ellipse_2;
    wg_2_x0(i) = xc_ellipse_2...
        - b0_ellipse_2 * sqrt(1 - (wg_2_z0(i) - zc_ellipse_2)^2 / a0_ellipse_2^2);
    wg_2_z1(i) = alpha * (zc_ellipse_2 - a1_ellipse_2)...
        + (1 - alpha) * zc_ellipse_2;
    wg_2_x1(i) = xc_ellipse_2...
        - b1_ellipse_2 * sqrt(1 - (wg_2_z1(i) - zc_ellipse_2)^2 / a1_ellipse_2^2);
end

imagesc(Eps_yy');
axis equal tight;

save 'Eps_yy.mat' Eps_yy wg_1_z0 wg_1_x0 wg_1_z1 wg_1_x1 ...
    wg_2_z0 wg_2_x0 wg_2_z1 wg_2_x1;
