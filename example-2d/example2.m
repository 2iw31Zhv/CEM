%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D FDTD Example (Ez Mode)
% Author: Ziwei Zhu
% Reference: http://emlab.utep.edu/ee5390fdtd.htm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear variables;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulation parameters
t_total = 2e-7; % s
x0 = 0.0; % m
x1 = 10.0; % m
y0 = 0.0; % m
y1 = 20.0; % m

% physical constant
c0 = 299792458; % m/s
epsilon0 = 8.854188e-12; % F/m

% device parameters
refractive_index_max = 1.0;
refractive_index_bc = 1.0;
refractive_index_source = 1.0;

devlen_min = 10.0;
grid_per_devlen = 4;

% EM paramters
freq_max = 5e8; % Hz
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

Nt = ceil( t_total ./ Dt);
Dt = t_total ./ Nt;

% 2x grid technique
Nx2 = 2 * Nx;
Ny2 = 2 * Ny;

Eps2 = ones(Nx2, Ny2);
Mu2 = ones(Nx2, Ny2);

% set device parameters here
Eps_device = 2.0;
Mu_device = 6.0;
device_center_x = 5.0;
device_center_y = 10.0;
device_radius = 3.0;

for nx = 1 : Nx2
    for ny = 1 : Ny2
        real_x = nx * res_x * 0.5;
        real_y = ny * res_y * 0.5;
        
        dist = sqrt((real_x - device_center_x)^2 ...
            + (real_y - device_center_y)^2);
        if dist <= device_radius
            Eps2(nx, ny) = Eps_device;
            Mu2(nx, ny) = Mu_device;
        end
    end
end


% extract the corresponding elements
Eps_zz = Eps2(1:2:Nx2, 1:2:Ny2);
Mu_xx = Mu2(1:2:Nx2, 2:2:Ny2);
Mu_yy = Mu2(2:2:Nx2, 1:2:Ny2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PML Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set PML parameters
N_layers_x0 = 20;
N_layers_x1 = 20;
N_layers_y0 = 20;
N_layers_y1 = 20;

Sigx2 = zeros(Nx2, Ny2);
% for nx = 1 : (2 * N_layers_x0)
%     pos_x = 2 * N_layers_x0 - nx + 1;
%     Sigx2(pos_x, :) = (0.5 * epsilon0 / Dt) * (nx / 2.0 / N_layers_x0)^3;
% end
% 
% for nx = 1 : (2 * N_layers_x1)
%     pos_x = Nx2 - 2 * N_layers_x1 + nx;
%     Sigx2(pos_x, :) = (0.5 * epsilon0 / Dt) * (nx / 2.0 / N_layers_x1)^3;
% end

Sigy2 = zeros(Nx2, Ny2);
for ny = 1 : (2 * N_layers_y0)
    pos_y = 2 * N_layers_y0 - ny + 1;
    Sigy2(:, pos_y) = (0.5 * epsilon0 / Dt) * (ny / 2.0 / N_layers_y0)^3;
end
for ny = 1 : (2 * N_layers_y1)
    pos_y = Ny2 - 2 * N_layers_y1 + ny;
    Sigy2(:, pos_y) = (0.5 * epsilon0 / Dt) * (ny / 2.0 / N_layers_y1)^3;
end

% now visualize it!
% subplot(1, 2, 1);
% imagesc(1:Nx2, 1:Ny2, Sigx2');
% subplot(1, 2, 2);
% imagesc(1:Nx2, 1:Ny2, Sigy2');

% compute PML update coeff
Sigx_Hx = Sigx2(1:2:Nx2, 2:2:Ny2);
Sigy_Hx = Sigy2(1:2:Nx2, 2:2:Ny2);

mHx0 = 1.0 / Dt + Sigy_Hx ./ (2.0 * epsilon0);
mHx1 = (1.0 / Dt - Sigy_Hx ./ (2.0 * epsilon0)) ./ mHx0;
mHx2 = (- c0 ./ Mu_xx) ./ mHx0;
mHx3 = (- c0 * Dt .* Sigx_Hx ./ (epsilon0 .* Mu_xx)) ./ mHx0;


Sigx_Hy = Sigx2(2:2:Nx2, 1:2:Ny2);
Sigy_Hy = Sigy2(2:2:Nx2, 1:2:Ny2);

mHy0 = 1.0 / Dt + Sigx_Hy ./ (2.0 * epsilon0);
mHy1 = (1.0 / Dt - Sigx_Hy ./ (2.0 * epsilon0)) ./ mHy0;
mHy2 = (-c0 ./ Mu_yy) ./ mHy0;
mHy3 = (-c0 * Dt .* Sigy_Hy ./ (epsilon0 .* Mu_yy)) ./ mHy0;


Sigx_Dz = Sigx2(1:2:Nx2, 1:2:Ny2);
Sigy_Dz = Sigy2(1:2:Nx2, 1:2:Ny2);

mDz0 = 1.0 / Dt + (Sigx_Dz + Sigy_Dz) ./ (2.0 * epsilon0)...
    + (Sigx_Dz .* Sigy_Dz) .* (Dt / 4.0 / epsilon0 / epsilon0);
mDz1 = (1.0 / Dt - (Sigx_Dz + Sigy_Dz) ./ (2.0 * epsilon0)...
    - (Sigx_Dz .* Sigy_Dz) .* (Dt / 4.0 / epsilon0 / epsilon0)) ./ mDz0;
mDz2 = c0 ./ mDz0;
mDz4 = -Dt / (epsilon0 * epsilon0) .* (Sigx_Dz .* Sigy_Dz) ./ mDz0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize HED field
Hx = zeros(Nx, Ny);
Hy = zeros(Nx, Ny);
Dz = zeros(Nx, Ny);
Ez = zeros(Nx, Ny);

CurlEx = zeros(Nx, Ny);
CurlEy = zeros(Nx, Ny);
CurlHz = zeros(Nx, Ny);

% initialize field integral
ICurlEx = zeros(Nx, Ny);
ICurlEy = zeros(Nx, Ny);
IDz = zeros(Nx, Ny);

% set source
tau = 0.5 ./ freq_max;
delay_time = 6.0 .* tau;
Generator = @(t)exp(-((t - delay_time) ./ tau).^2);

t_array = Dt * (1 : Nt);
source_array = Generator(t_array);
source_x = 5;
source_y = 1.0;

Ns_x = round((source_x - x0) ./ res_x);
Ns_y = round((source_y - y0) ./ res_y);


e_source = source_array;
h_source = sqrt( Eps_zz(Ns_x, Ns_y) / Mu_xx(Ns_x, Ns_y-1)) ...
    * Generator(t_array + 0.5 * Dt...
    + refractive_index_source * res_y / 2.0 / c0);

% for visualization
x_array = linspace(x0, x1, Nx);
y_array = linspace(y0, y1, Ny);

% Ez_recorder = zeros(Nt, 1);

% initialize movie recorder
movie_name = '2D_FDTD_Example';
delete(strcat(movie_name, '.mp4'));
vidObj = VideoWriter(movie_name, 'MPEG-4');
open(vidObj);

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
 % set device drawing
phi = linspace(0, 2*pi, 100);
profile_x = device_radius * cos(phi) + device_center_x;
profile_y = device_radius * sin(phi) + device_center_y;
device_color = [0.7 0.7 0.7];

% HDE Algorithm main loop
for T = 1 : Nt
    
    % evaluate Curl Ex
    for nx = 1 : Nx
        for ny = 1 : (Ny - 1)
            CurlEx(nx, ny) = (Ez(nx, ny+1) - Ez(nx, ny)) ./ res_y;
        end
        % handle boundary condition
        CurlEx(nx, Ny) = (0.0 - Ez(nx, Ny)) ./ res_y;
    end
    
    % TF/SF for E
    CurlEx(:, Ns_y-1) = CurlEx(:, Ns_y-1) - e_source(T) / res_y;
    
    ICurlEx = ICurlEx + CurlEx;
    
    % evaluate Curl Ey
    for ny = 1 : Ny
        for nx = 1 : (Nx - 1)
            CurlEy(nx, ny) = -(Ez(nx+1, ny) - Ez(nx, ny)) ./ res_x;
        end
        % handle boundary condition
        % CurlEy(Nx, ny) = -(0.0 - Ez(Nx, ny)) ./ res_x;
        CurlEy(Nx, ny) = -(Ez(1, ny) - Ez(nx, ny)) ./ res_x;
    end
    
    ICurlEy = ICurlEy + CurlEy;
    
    % Hx, Hy <- Curl Ex, Curl Ey
    Hx = mHx1 .* Hx + mHx2 .* CurlEx + mHx3 .* ICurlEx;
    Hy = mHy1 .* Hy + mHy2 .* CurlEy + mHy3 .* ICurlEy;
    
    % evaluate Curl Hz
    CurlHz(1, 1) = (Hy(1, 1) - Hy(Nx, 1)) ./ res_x...
                - (Hx(1, 1) - 0.0) ./ res_y;
    for ny = 2 : Ny
        CurlHz(1, ny) = (Hy(1, ny) - Hy(Nx, ny)) ./ res_x...
                - (Hx(1, ny) - Hx(1, ny-1)) ./ res_y;
    end
    
    for nx = 2 : Nx
        CurlHz(nx, 1) = (Hy(nx, 1) - Hy(nx-1, 1)) ./ res_x...
                - (Hx(nx, 1) - 0.0) ./ res_y;
        for ny = 2 : Ny
            CurlHz(nx, ny) = (Hy(nx, ny) - Hy(nx-1, ny)) ./ res_x...
                - (Hx(nx, ny) - Hx(nx, ny-1)) ./ res_y;
        end
    end
    
    % TF/SF for Hx
    CurlHz(:, Ns_y) = CurlHz(:, Ns_y) + h_source(T) / res_y;
    
    % Dz <- Curl Hz
    % Dz = Dz + c0 .* Dt .* CurlHz;
    Dz = mDz1 .* Dz + mDz2 .* CurlHz + mDz4 .* IDz;
    
    IDz = IDz + Dz;
    
    % simple soft source (not recommended)
    % Dz(Ns_x, Ns_y) = Dz(Ns_x, Ns_y) + source_array(T);
    
    % Ez <- Dz
    Ez = Dz ./ Eps_zz;
    
    % visualize the results
    if mod(T, 5) == 0
        shift = 8;
        clf;
        
        subplot(1, 3, 1);
        fill(profile_x, profile_y, device_color);
        hold on;
        imagesc(x_array, y_array, field_log_normalize(Ez, shift)');
        caxis([-shift, shift]);
        axis equal tight;
        title('Ez');
        alpha(0.3);
        
        subplot(1, 3, 2);
        fill(profile_x, profile_y, device_color);
        hold on;
        imagesc(x_array, y_array, field_log_normalize(Hx, shift)');
        caxis([-shift, shift]);
        axis equal tight;
        title('Hx');
        alpha(0.3);
        
        
        subplot(1, 3, 3);
        fill(profile_x, profile_y, device_color);
        hold on;
        imagesc(x_array, y_array, field_log_normalize(Hy, shift)');
        caxis([-shift, shift]);
        axis equal tight;
        title('Hy');
        alpha(0.3);
                
        t = Dt .*T;
        title_str = sprintf('2D FDTD Example (Ez Mode), t = %.3e s', t);
        subtitle(title_str);
        
        F = getframe(gcf);
        writeVideo(vidObj, F);
        % pause(0.001);
    end
end

close(vidObj);