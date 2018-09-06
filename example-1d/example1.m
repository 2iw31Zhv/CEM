%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Example code for 1D FDTD
% ===============================================
% Author: Ziwei Zhu
% Reference: http://emlab.utep.edu/ee5390fdtd.htm
% Units: GI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

% physical constant
c0 = 299792458; % m/s
mu = 1; % H/m
epsilon = 1; % F/m

% initalize physical domain
length_z = 40.0; % m


% compute the resolution of the grid
% wavelength part
N_wave = 20; 
max_freq = 5e7; % Hz
n_max = 1.0;
min_wavelen = c0 ./ max_freq ./ n_max;
res_wave = min_wavelen ./ N_wave;

% device part
N_device = 4;
min_feature = 0.1; % m
res_device = min_feature ./ N_device;

res_z = min(res_wave, res_device);

Nz = ceil(length_z / res_z);
res_z = length_z / Nz;

z_array_Ey = res_z .* (0 : (Nz-1));
z_array_Hx = z_array_Ey + 0.5 .* res_z;

% compute the source pulse
tau = 0.5 / max_freq;
delay_time = 6 * tau;


% evaluate time step
n_bc = 1.0;
Dt = n_bc .* res_z ./ 2.0 ./ c0;


% evaluate total time
t_prop = n_max .* length_z ./ c0;
t_total = 12 * tau + 5 * t_prop;

n_steps = ceil(t_total ./ Dt);

% initialize Ey / Hx field
Hx = zeros(Nz, 1);
Ey = zeros(Nz, 1);
E1 = 0;
E2 = 0;
H1 = 0;
H2 = 0;


% initialize update coefficients
m_Ey = c0 .* Dt ./ epsilon * ones(Nz, 1);
m_Hx = c0 .* Dt ./ mu * ones(Nz, 1);

% evaluate source
t_array = Dt * (1 : n_steps);
e_source = exp(-((t_array - delay_time) ./ tau).^2);
pos_source = round(0.4 * Nz);

% compute source correction term
e_correction = e_source;
n_source = 1.0;
t_shift = 0.5 * Dt - n_source * res_z / 2.0 / c0;
h_correction = -sqrt(epsilon / mu) ...
    * exp(-((t_array - delay_time + t_shift) ./ tau).^2);

% initialize movie recorder
movie_name = '1D_FDTD_Example2';
delete(strcat(movie_name, '.mp4'));
vidObj = VideoWriter(movie_name, 'MPEG-4');
open(vidObj);

% update the electric and magnetic field
for n = 1 : n_steps
    
    % record before H
    H2 = H1;
    H1 = Hx(1);
    
    % update H <- E
    for nz = 1 : (Nz-1)
       Hx(nz) = Hx(nz) + m_Hx(nz) .* (Ey(nz+1) - Ey(nz)) ./ res_z;
    end
    
    % H source correction
    Hx(pos_source - 1) = Hx(pos_source - 1)...
        -  m_Hx(pos_source - 1) .* e_correction(n) ./ res_z;
    
    % handle perfect boundary
    Hx(Nz) = Hx(Nz) + m_Hx(Nz) .* (E2 - Ey(Nz)) ./ res_z;
    % or
    % pure conductor boundary
    % Hx(Nz) = Hx(Nz) + m_Hx(Nz) .* (0 - Ey(Nz)) ./ res_z;
    
    % record before E
    E2 = E1;
    E1 = Ey(Nz);
    
    % handle perfect boundary
    Ey(1) = Ey(1) + m_Ey(1) .* (Hx(1) - H2) ./ res_z;
    % or
    % pure conductor boundary
    % Ey(1) = Ey(1) + m_Ey(1) .* (Hx(1) - 0) ./ res_z;
    
    % update E <- H
    for nz = 2 : Nz
       Ey(nz) = Ey(nz) + m_Ey(nz) .* (Hx(nz) - Hx(nz-1)) ./ res_z;
    end
    
    % E source correction
    Ey(pos_source) = Ey(pos_source) - m_Ey(pos_source)...
        * h_correction(n) ./ res_z;
    
    if mod(n, 50) == 1
        % save the frame to the movie
        clf;
        
        plot(z_array_Ey, Ey, '-b', 'LineWidth', 3);
        hold on;
        plot(z_array_Hx, Hx, '-r', 'LineWidth', 2);
        axis([0, length_z, -2, 2]);
        legend('E_y', 'H_x');
        xlabel('position (m)');
        ylabel('Intensity (V/m)');
        t = Dt * n;
        title_str = sprintf('1D FDTD Example, t = %.3e s', t);
        title(title_str);
        
        F = getframe(gcf);
        writeVideo(vidObj, F);
        
    end
end

close(vidObj);