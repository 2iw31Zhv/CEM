%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D FDTD Example (Ez Mode)
% Author: Ziwei Zhu
% Reference: http://emlab.utep.edu/ee5390fdtd.htm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulation parameters
t_total = 10.0; % s
x0 = 0.0; % m
x1 = 10.0; % m
y0 = 0.0; % m
y1 = 10.0; % m

% physical constant
c0 = 299792458; % m/s

% device parameters
refractive_index_max = 1.0;
refractive_index_bc = 1.0;

devlen_min = 10.0;
grid_per_devlen = 4;

% EM paramters
freq_max = 5e7; % Hz
grid_per_wavelen = 20;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evaluate resolution
wavelen_min = c0 ./ freq_max ./ refractive_index_max;
res_wave = wavelen_min ./ grid_per_wavelen;
res_dev = devlen_min ./ grid_per_devlen;
res = min(res_wave, res_dev);

% evaluate grid number
Nx = ceil((x1 - x0) ./ res);
Ny = ceil((y1 - y0) ./ res);
res_x = (x1 - x0) ./ Nx;
res_y = (y1 - y0) ./ Ny;

% evaluate time step
res_min = min(res_x, res_y);
Dt = refractive_index_bc .* res_min ./ 2.0 ./ c0;

% initialize HED field
