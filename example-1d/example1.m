%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Example code for 1D FDTD
% ===============================================
% Author: Ziwei Zhu
% Reference: http://emlab.utep.edu/ee5390fdtd.htm
% Units: GI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% physical constant
c0 = 299792458; % m/s

% initalize physical domain
length_z = 10.0; % m


% compute the resolution of the grid
% wavelength part
N_wave = 20; 
max_freq = 5e14; % Hz
max_refractive_index = 3.0;
min_wavelen = c0 ./ max_freq ./ max_refractive_index;
res_wave = min_wavelen ./ N_wave;

% device part
N_device = 4;
min_feature = 0.1; % m
res_device = min_feature ./ N_device;

res_z = min(res_wave, res_device);

Nz = ceil(length_z / res_z);
res_z = length_z / Nz;

% compute the source pulse
tau = 0.5 / max_freq;
delay_time = 6 * tau;


% initialize Ey / Hx field
Hx = zeros(Nz, 1);
Ey = zeros(Nz, 1);


