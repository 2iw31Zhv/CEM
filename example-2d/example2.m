%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D FDTD Example (Ez Mode)
% Author: Ziwei Zhu
% Reference: http://emlab.utep.edu/ee5390fdtd.htm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulation parameters
t_total = 1e-6; % s
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

% set permitivity and permeability
Mu_xx = ones(Nx, Ny);
Mu_yy = ones(Nx, Ny);
Eps_zz = ones(Nx, Ny);

% initialize HED field
Hx = zeros(Nx, Ny);
Hy = zeros(Nx, Ny);
Dz = zeros(Nx, Ny);
Ez = zeros(Nx, Ny);

CurlEx = zeros(Nx, Ny);
CurlEy = zeros(Nx, Ny);
CurlHz = zeros(Nx, Ny);

% set source
tau = 0.5 ./ freq_max;
delay_time = 6.0 .* tau;
Generator = @(t)exp(-((t - delay_time) ./ tau).^2);

source_x = 2.5;
source_y = 2.5;

Ns_x = round((source_x - x0) ./ res_x);
Ns_y = round((source_y - y0) ./ res_y);


% for visualization
x_array = linspace(x0, x1, Nx);
y_array = linspace(y0, y1, Ny);

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
    
    % evaluate Curl Ey
    for ny = 1 : Ny
        for nx = 1 : (Nx - 1)
            CurlEy(nx, ny) = -(Ez(nx+1, ny) - Ez(nx, ny)) ./ res_x;
        end
        % handle boundary condition
        CurlEy(Nx, ny) = -(0.0 - Ez(Nx, ny)) ./ res_x;
    end
    
    % Hx, Hy <- Curl Ex, Curl Ey
    Hx = Hx - (c0 .* Dt ./ Mu_xx) .* CurlEx;
    Hy = Hy - (c0 .* Dt ./ Mu_yy) .* CurlEy;
    
    % evaluate Curl Hz
    CurlHz(1, 1) = (Hy(1, 1) - 0.0) ./ res_x...
                - (Hx(1, 1) - 0.0) ./ res_y;
    for ny = 2 : Ny
        CurlHz(1, ny) = (Hy(1, ny) - 0.0) ./ res_x...
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
    
    
    % Dz <- Curl Hz
    Dz = Dz + c0 .* Dt .* CurlHz;
    
    % simple soft source
    t_now = Dt .* T;
    Dz(Ns_x, Ns_y) = Dz(Ns_x, Ns_y) + Generator(t_now);
    
    % Ez <- Dz
    Ez = Dz ./ Eps_zz;
    
    % visualize the results
    if mod(T, 1) == 0
        
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
        
        imagesc(x_array, y_array, Ez');
        t = Dt .*T;
        title_str = sprintf('2D FDTD Example (Ez Mode), t = %.3e s', t);
        title(title_str);
        pause(0.01);
    end
end