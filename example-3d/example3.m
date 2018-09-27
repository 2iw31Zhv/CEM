%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D FDTD Example
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
x0 = 0.0;
x1 = 10.0;
y0 = 0.0;
y1 = 10.0;
z0 = 0.0;
z1 = 40.0;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main HDE Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

