%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D FDTD Example
% =========================================================================
% Description:
% Assume periodic boundary condition for both x and y directions
% and perfect match layer along the z direction
% =========================================================================
% Author: Ziwei Zhu
% Reference: http://emlab.utep.edu/ee5390fdtd.htm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear variables;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% physical constant
c0 = 299792458; % m/s
epsilon0 = 8.854188e-12; % F/m

% simulation size
x0 = -5.0;
x1 = 5.0;
y0 = -5.0;
y1 = 5.0;
z0 = -20.0;
z1 = 20.0;

t_total = 2e-7;

% device parameters
devlen_min = 5.0;
grid_per_devlen = 4;


% EM parameters
freq_max = 5e7;
grid_per_wavelen = 20;

% required dependent parameters
refractive_index_max = 1.0;
refractive_index_bc = 1.0;
refractive_index_source = 1.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evaluate resolution
wavelen_min = c0 ./ freq_max ./ refractive_index_max;
res_wave = wavelen_min ./ grid_per_wavelen;
res_dev = devlen_min ./ grid_per_devlen;
res = min(res_wave, res_dev);

% evaluate grid number
Nx = ceil((x1 - x0) ./ res);
Ny = ceil((y1 - y0) ./ res);
Nz = ceil((z1 - z0) ./ res);

res_x = (x1 - x0) ./ Nx;
res_y = (y1 - y0) ./ Ny;
res_z = (z1 - z0) ./ Nz;

res_min = min([res_x res_y res_z]);


% evaluate time step
Dt = refractive_index_bc .* res_min ./ 2.0 ./ c0;

% evaluate iter number
Nt = ceil( t_total ./ Dt);
Dt = t_total ./ Nt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Large Scale Parameter Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2x grid technique
Nx2 = 2 * Nx;
Ny2 = 2 * Ny;
Nz2 = 2 * Nz;

Eps2 = ones(Nx2, Ny2, Nz2);
Mu2 = ones(Nx2, Ny2, Nz2);

% TODO 
% set device epsilon and mu here!

% extract the corresponding elements
Eps_xx = Eps2(2:2:Nx2, 1:2:Ny2, 1:2:Nz2);
Eps_yy = Eps2(1:2:Nx2, 2:2:Ny2, 1:2:Nz2);
Eps_zz = Eps2(1:2:Nx2, 1:2:Ny2, 2:2:Nz2);

Mu_xx = Mu2(1:2:Nx2, 2:2:Ny2, 2:2:Nz2);
Mu_yy = Mu2(2:2:Nx2, 1:2:Ny2, 2:2:Nz2);
Mu_zz = Mu2(2:2:Nx2, 2:2:Ny2, 1:2:Nz2);


% set PML parameters
N_layers_z0 = 20;
N_layers_z1 = 20;

Sigz2 = zeros(Nx2, Ny2, Nz2);
for nz = 1 : (2 * N_layers_z0)
    pos_z = 2 * N_layers_z0 - nz + 1;
    Sigz2(:, :, pos_z) = (0.5 * epsilon0 / Dt) * (nz / 2.0 / N_layers_z0)^3;
end
for nz = 1 : (2 * N_layers_z1)
    pos_z = Nz2 - 2 * N_layers_z1 + nz;
    Sigz2(:, :, pos_z) = (0.5 * epsilon0 / Dt) * (nz / 2.0 / N_layers_z1)^3;
end

x_array = x0 + res_x * (1 : Nx);
y_array = y0 + res_y * (1 : Ny);
z_array = z0 + res_z * (1 : Nz);
[Y, X, Z] = meshgrid(y_array, x_array, z_array);

Sigz_Hx = Sigz2(1:2:Nx2, 2:2:Ny2, 2:2:Nz2);
Sigz_Hy = Sigz2(2:2:Nx2, 1:2:Ny2, 2:2:Nz2);
Sigz_Hz = Sigz2(2:2:Nx2, 2:2:Ny2, 1:2:Nz2);

Sigz_Dx = Sigz2(2:2:Nx2, 1:2:Ny2, 1:2:Nz2);
Sigz_Dy = Sigz2(1:2:Nx2, 2:2:Ny2, 1:2:Nz2);
Sigz_Dz = Sigz2(1:2:Nx2, 1:2:Ny2, 2:2:Nz2);

% evaluate update coeff
mHx0 = 1.0 / Dt + Sigz_Hx ./ (2.0 * epsilon0);
mHx1 = (1.0 / Dt - Sigz_Hx ./ (2.0 * epsilon0)) ./ mHx0;
mHx2 = (- c0 ./ Mu_xx) ./ mHx0;

mHy0 = 1.0 / Dt + Sigz_Hy ./ (2.0 * epsilon0);
mHy1 = (1.0 / Dt - Sigz_Hy ./ (2.0 * epsilon0)) ./ mHy0;
mHy2 = (-c0 ./ Mu_yy) ./ mHy0;

mHz2 = - c0 * Dt ./ Mu_zz;
mHz3 = - c0 * Dt * Dt / epsilon0 * Sigz_Hz ./ Mu_zz;

mDx0 = 1.0 / Dt + Sigz_Dx ./ (2.0 * epsilon0);
mDx1 = (1.0 / Dt - Sigz_Dx ./ (2.0 * epsilon0)) ./ mDx0;
mDx2 = c0 ./ mDx0;

mDy0 = 1.0 / Dt + Sigz_Dy ./ (2.0 * epsilon0);
mDy1 = (1.0 / Dt - Sigz_Dy ./ (2.0 * epsilon0)) ./ mDy0;
mDy2 = c0 ./ mDy0;

mDz2 = c0 * Dt * ones(Nx, Ny, Nz);
mDz3 = c0 * Dt * Dt / epsilon0 * Sigz_Dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main HDE Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize field, curl and integrator
Hx = zeros(Nx, Ny, Nz);
Hy = zeros(Nx, Ny, Nz);
Hz = zeros(Nx, Ny, Nz);

Dx = zeros(Nx, Ny, Nz);
Dy = zeros(Nx, Ny, Nz);
Dz = zeros(Nx, Ny, Nz);

Ex = zeros(Nx, Ny, Nz);
Ey = zeros(Nx, Ny, Nz);
Ez = zeros(Nx, Ny, Nz);

CurlEx = zeros(Nx, Ny, Nz);
CurlEy = zeros(Nx, Ny, Nz);
CurlEz = zeros(Nx, Ny, Nz);

ICurlEz = zeros(Nx, Ny, Nz);

CurlHx = zeros(Nx, Ny, Nz);
CurlHy = zeros(Nx, Ny, Nz);
CurlHz = zeros(Nx, Ny, Nz);

ICurlHz = zeros(Nx, Ny, Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the colormap from http://emlab.utep.edu/ee5390fdtd.htm
CMAP = zeros(256,3);
c1 = [0 0 1]; %blue
c2 = [1 1 1]; %white
c3 = [1 0 0]; %red
for nc = 1 : 128
    f = (nc - 1)/128;
    c = (1 - sqrt(f))*c1 + sqrt(f)*c2;
    CMAP(nc,:) = c;
    c = (1 - f^2)*c2 + f^2*c3;
    CMAP(128+nc,:) = c;
end
colormap(CMAP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% main loop
for T = 1 : Nt
    % =====================================================================
    % here we use circshift(..., -1) to get the next element according to 
    % the periodic boundary condition
    % =====================================================================
    % evaluate Curl Ex
    CurlEx = (circshift(Ez, [0, -1, 0]) - Ez) ./ res_y...
        - (circshift(Ey, [0, 0, -1]) - Ey) ./ res_z;
    
    % evaluate Curl Ey
    CurlEy = (circshift(Ex, [0, 0, -1]) - Ex) ./ res_z...
        - (circshift(Ez, [-1, 0, 0]) - Ez) ./ res_x;
    
    % evaluate Curl Ez
    CurlEz = (circshift(Ey, [-1, 0, 0]) - Ey) ./ res_x...
        - (circshift(Ex, [0, -1, 0]) - Ex) ./ res_y;
    
    ICurlEz = ICurlEz + CurlEz;
    
    % evaluate Curl Hx
    CurlHx = (Hz - circshift(Hz, [0, 1, 0])) ./ res_y...
        - (Hy - circshift(Hy, [0, 0, 1])) ./ res_z;
    
    % evaluate Curl Hy
    CurlHy = (Hx - circshift(Hx, [0, 0, 1])) ./ res_z...
        - (Hz - circshift(Hz, [1, 0, 0])) ./ res_x;
    
    % evaluate Curl Hz
    CurlHz = (Hy - circshift(Hy, [1, 0, 0])) ./ res_x...
        - (Hx - circshift(Hx, [0, 1, 0])) ./ res_y;
    
    ICurlHz = ICurlHz + CurlHz;
    
    % update Hx, Hy, Hz
    Hx = mHx1 .* Hx + mHx2 .* CurlEx;
    Hy = mHy1 .* Hy + mHy2 .* CurlEy;
    Hz = Hz + mHz2 .* CurlEz + mHz3 .* ICurlEz;
    
    % update Dx, Dy, Dz
    Dx = mDx1 .* Dx + mDx2 .* CurlHx;
    Dy = mDy1 .* Dy + mDy2 .* CurlHy;
    Dz = Dz + mDz2 .* CurlHz + mDz3 .* ICurlHz;
    
    % update Ex, Ey, Ez
    Ex = Dx ./ Eps_xx;
    Ey = Dy ./ Eps_yy;
    Ez = Dz ./ Eps_zz;
    
    % visualize
    if mod(T, 20) == 0
        slice(Y, X, -Z, Ex, 0, 0, 0);
        axis equal tight off;
        shading interp;
        grid on;
        pause(0.01);
    end
end


