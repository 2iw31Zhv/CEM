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
t_total = 1e-6; % s
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
if mod(Nx, 2) == 0
    Nx = Nx + 1;
end

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
% figure('units','normalized','outerposition',[0 0 1 1]);
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

% compute steady state Ez field using fourier transform
N_record_ref = Ns_y - 2;
N_record_trn = Ny - N_layers_y1 - 1;

n_freq = 100;
freq_array = linspace(0, freq_max, n_freq);
freq_kernel = exp(-1i * 2 * pi * freq_array * Dt);

REF = zeros(Nx, n_freq);
TRN = zeros(Nx, n_freq);
SRC = zeros(Nx, n_freq);

DE_r_allmodes = zeros(n_freq, 1);
DE_t_allmodes = zeros(n_freq, 1);

Eref = zeros(Nx, n_freq);
Etrn = zeros(Nx, n_freq);
    
% HDE Algorithm main loop
for T = 1 : Nt
    % evaluate Curl Ex
    CurlEx = (circshift(Ez, [0 -1]) - Ez) / res_y;
    CurlEx(:, Ny) = (0.0 - Ez(:, Ny)) ./ res_y;
    
    % TF/SF for E
    CurlEx(:, Ns_y-1) = CurlEx(:, Ns_y-1) - e_source(T) / res_y;
    
    ICurlEx = ICurlEx + CurlEx;
    
    % evaluate Curl Ey 
    CurlEy = - (circshift(Ez, [-1 0]) - Ez) / res_x;
    
    ICurlEy = ICurlEy + CurlEy;
    
    % Hx, Hy <- Curl Ex, Curl Ey
    Hx = mHx1 .* Hx + mHx2 .* CurlEx + mHx3 .* ICurlEx;
    Hy = mHy1 .* Hy + mHy2 .* CurlEy + mHy3 .* ICurlEy;
    
    % evaluate Curl Hz
    CurlHz = (Hy - circshift(Hy, [1 0])) / res_x...
        - (Hx - circshift(Hx, [0 1])) / res_y;
    CurlHz(:, 1) = (Hy(:, 1) - circshift(Hy(:, 1), [1 0])) / res_x...
        - (Hx(:, 1) - 0.0) / res_y;
    
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
    
    for nf = 1 : n_freq
        REF(:, nf) = REF(:, nf) + Dt * (freq_kernel(nf).^T) * Ez(:, N_record_ref);
        TRN(:, nf) = TRN(:, nf) + Dt * (freq_kernel(nf).^T) * Ez(:, N_record_trn);
        SRC(:, nf) = SRC(:, nf) + Dt * (freq_kernel(nf).^T) * e_source(T);
    end
    
    % get reflectance and transmittance for different Bloch modes
    n_inc = 1.0;
    n_record_ref = 1.0;
    n_record_trn = 1.0;

    for nfreq = 1 : n_freq
        
        k0 = 2 * pi * freq_array(nfreq) / c0;
        % for normal incident light
        kyinc = n_inc * k0;

        m = (- floor(Nx/2) : floor(Nx/2))';
        kxm = -2 * pi * m / (x1 - x0);

        ky_ref = sqrt((k0 * n_record_ref)^2 - kxm.^2);
        ky_trn = sqrt((k0 * n_record_trn)^2 - kxm.^2);

        Eref(:, nfreq) = REF(:, nfreq) ./ SRC(:, nfreq);
        Etrn(:, nfreq) = TRN(:, nfreq) ./ SRC(:, nfreq);

        Eref(:, nfreq) = fftshift(fft(Eref(:, nfreq))) / Nx;
        Etrn(:, nfreq) = fftshift(fft(Etrn(:, nfreq))) / Nx;

        DE_ref = abs(Eref).^2 .* real(ky_ref / kyinc);
        DE_trn = abs(Etrn).^2 .* real(ky_trn / kyinc);

        DE_r_allmodes(nfreq) = sum(DE_ref(:, nfreq));
        DE_t_allmodes(nfreq) = sum(DE_trn(:, nfreq));
    end

    DE_s_allmodes = DE_r_allmodes + DE_t_allmodes;

    % visualize the results
    if mod(T, 10) == 0
        shift = 8;
        clf;
        
        subplot(2, 3, 1);
        fill(profile_x, profile_y, device_color);
        hold on;
        imagesc(x_array, y_array, field_log_normalize(Ez, shift)');
        caxis([-shift, shift]);
        axis equal tight;
        xlabel('Ez');
        %title('Ez');
        alpha(0.3);
        
        subplot(2, 3, 2);
        fill(profile_x, profile_y, device_color);
        hold on;
        imagesc(x_array, y_array, field_log_normalize(Hx, shift)');
        caxis([-shift, shift]);
        axis equal tight;
        xlabel('Hx');
        %title('Hx');
        alpha(0.3);
        
        subplot(2, 3, 3);
        fill(profile_x, profile_y, device_color);
        hold on;
        imagesc(x_array, y_array, field_log_normalize(Hy, shift)');
        caxis([-shift, shift]);
        axis equal tight;
        xlabel('Hy');
        %title('Hy');
        alpha(0.3);
        
        subplot(2, 3, 4:6);
        plot(freq_array, DE_r_allmodes, '-r', 'LineWidth', 3);
        hold on;
        plot(freq_array, DE_t_allmodes, '-g', 'LineWidth', 2);
        hold on;
        plot(freq_array, DE_s_allmodes, '-k', 'LineWidth', 1);
        axis([freq_array(2), freq_max, 0, 1.5]);
        xlabel('Frequency (Hz)');
        legend('reflectance', 'transmittance', 'total');
              
        t = Dt .* T;
        title_str = sprintf('2D FDTD Example (Ez Mode), t = %.3e s', t);
        subtitle(title_str);
        
        F = getframe(gcf);
        writeVideo(vidObj, F);
        % pause(0.001);
    end
end
close(vidObj);

save 'transmittance.mat' DE_trn
