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

t_total = 2e-6;

% device parameters
devlen_min = 4.0;
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
wavelen_min = c0 / freq_max / refractive_index_max;
res_wave = wavelen_min / grid_per_wavelen;
res_dev = devlen_min / grid_per_devlen;
res = min(res_wave, res_dev);

% evaluate grid number
Nx = ceil((x1 - x0) / res);
Ny = ceil((y1 - y0) / res);
Nz = ceil((z1 - z0) / res);

res_x = (x1 - x0) / Nx;
res_y = (y1 - y0) / Ny;
res_z = (z1 - z0) / Nz;

res_min = min([res_x res_y res_z]);

% evaluate time step
Dt = refractive_index_bc * res_min / 2.0 / c0;

% refine time step according to source
tau = 0.5 / freq_max;
N_source = 20;
res_source = tau / N_source;
Dt = min(Dt, res_source);

% evaluate iter number
Nt = ceil( t_total / Dt);
Dt = t_total / Nt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Large Scale Parameter Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2x grid technique
Nx2 = 2 * Nx;
Ny2 = 2 * Ny;
Nz2 = 2 * Nz;

Eps2 = ones(Nx2, Ny2, Nz2);
Mu2 = ones(Nx2, Ny2, Nz2);

% set device epsilon and mu here!

dev_x0 = -2.0;
dev_x1 = 2.0;
dev_y0 = -3.0;
dev_y1 = 3.0;
dev_z0 = -4.0;
dev_z1 = 4.0;

dev_eps = 2.0;
dev_mu = 6.0;

for nx = 1 : Nx2
    for ny = 1 : Ny2
        for nz = 1 : Nz2
            x = x0 + 0.5 * res_x * nx;
            y = y0 + 0.5 * res_y * ny;
            z = z0 + 0.5 * res_z * nz;
            
            if (dev_x0 <= x && x <= dev_x1...
                    && dev_y0 <= y && y <= dev_y1...
                    && dev_z0 <= z && z <= dev_z1)
                Eps2(nx, ny, nz) = dev_eps;
                Mu2(nx, ny, nz) = dev_mu;
            end
        end
    end
end

L = dev_x1 - dev_x0;
W = dev_y1 - dev_y0;
H = dev_z1 - dev_z0;

A = [dev_x0;dev_y0;dev_z0];
B = A + [L;0;0];
C = B + [0;W;0];
D = A + [0;W;0];

r1 = repmat(A,1,5);
r2 = [A,B,C,D,A];
r3 = r2 + repmat([0;0;H],1,5);
r4 = repmat(r3(:,1),1,5);
devP=[r1;r2;r3;r4];

devX = devP(1:3:end,:);
devY = devP(2:3:end,:);
devZ = devP(3:3:end,:);

dev_x = devX(2:3,:);
dev_y = devY(2:3,:);
dev_z = devZ(2:3,:);

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
mHx2 = (-c0 ./ Mu_xx) ./ mHx0;

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

mDz2 = c0 * Dt;
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

% set source
delay_time = 6.0 * tau;
Generator = @(t)exp(-((t - delay_time) ./ tau).^2);

t_array = Dt * (1 : Nt);
e_source_array = Generator(t_array);
h_source_array = Generator(t_array + 0.5 * Dt ...
    + refractive_index_source * res_z / 2.0 / c0);

source_x = 0.0;
source_y = 0.0;
source_z = -10.0;

source_angle = pi / 4;
Px = cos(source_angle);
Py = sin(source_angle);

Ns_x = round((source_x - x0) ./ res_x);
Ns_y = round((source_y - y0) ./ res_y);
Ns_z = round((source_z - z0) ./ res_z);

Ex_source = Px * e_source_array;
Ey_source = Py * e_source_array;

impedance_ref_x = sqrt(Mu_xx(Ns_x, Ns_y, Ns_z) / Eps_xx(Ns_x, Ns_y, Ns_z));
impedance_ref_y = sqrt(Mu_yy(Ns_x, Ns_y, Ns_z) / Eps_yy(Ns_x, Ns_y, Ns_z));

Hx_source = -Py * h_source_array / impedance_ref_x;
Hy_source = Px * h_source_array / impedance_ref_y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the colormap from http://emlab.utep.edu/ee5390fdtd.htm
figure('units','normalized','outerposition',[0 0 1 1]);
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

record = zeros([Nt, 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize movie recorder
movie_name = '3D_FDTD_Example';
delete(strcat(movie_name, '.mp4'));
vidObj = VideoWriter(movie_name, 'MPEG-4');
open(vidObj);

% main loop
for T = 1 : Nt
    % =====================================================================
    % here we use circshift(..., -1) to get the next element according to 
    % the periodic boundary condition
    % =====================================================================
    % evaluate Curl Ex
    CurlEx = (circshift(Ez, [0, -1, 0]) - Ez) ./ res_y...
        - (circshift(Ey, [0, 0, -1]) - Ey) ./ res_z;
    % handle z boundary
    CurlEx(:, :, Nz) = (circshift(Ez(:, :, Nz), [0, -1]) - Ez(:, :, Nz)) ./ res_y...
        - (0.0 - Ey(:, :, Nz)) ./ res_z;
    
    CurlEx(:, :, Ns_z-1) = CurlEx(:, :, Ns_z-1)...
        + Ey_source(T) / res_z;
    
    % evaluate Curl Ey
    CurlEy = (circshift(Ex, [0, 0, -1]) - Ex) ./ res_z...
        - (circshift(Ez, [-1, 0, 0]) - Ez) ./ res_x;
    % handle z boundary
    CurlEy(:, :, Nz) = (0.0 - Ex(:, :, Nz)) ./ res_z...
        - (circshift(Ez(:, :, Nz), [-1, 0]) - Ez(:, :, Nz)) ./ res_x;
    
    CurlEy(:, :, Ns_z-1) = CurlEy(:, :, Ns_z-1)...
        - Ex_source(T) / res_z;
    
    % evaluate Curl Ez
    CurlEz = (circshift(Ey, [-1, 0, 0]) - Ey) ./ res_x...
        - (circshift(Ex, [0, -1, 0]) - Ex) ./ res_y;
    
    ICurlEz = ICurlEz + CurlEz;
    
        
    % update Hx, Hy, Hz
    Hx = mHx1 .* Hx + mHx2 .* CurlEx;
    Hy = mHy1 .* Hy + mHy2 .* CurlEy;
    Hz = Hz + mHz2 .* CurlEz + mHz3 .* ICurlEz;
   
    % evaluate Curl Hx
    CurlHx = (Hz - circshift(Hz, [0, 1, 0])) ./ res_y...
        - (Hy - circshift(Hy, [0, 0, 1])) ./ res_z;
    % handle z boudary
    CurlHx(:, :, 1) = (Hz(:, :, 1) - circshift(Hz(:, :, 1), [0, 1])) ./ res_y...
        - (Hy(:, :, 1) - 0.0) ./ res_z;
    CurlHx(:, :, Ns_z) = CurlHx(:, :, Ns_z)...
        + Hy_source(T) / res_z;
    
    % evaluate Curl Hy
    CurlHy = (Hx - circshift(Hx, [0, 0, 1])) ./ res_z...
        - (Hz - circshift(Hz, [1, 0, 0])) ./ res_x;
    % handle z boundary
    CurlHy(:, :, 1) = (Hx(:, :, 1) - 0.0) ./ res_z...
        - (Hz(:, :, 1) - circshift(Hz(:, :, 1), [1, 0])) ./ res_x;
    
    CurlHy(:, :, Ns_z) = CurlHy(:, :, Ns_z)...
        - Hx_source(T) / res_z;
    
    
    ICurlHz = ICurlHz + CurlHz;
    
    % evaluate Curl Hz
    CurlHz = (Hy - circshift(Hy, [1, 0, 0])) ./ res_x...
        - (Hx - circshift(Hx, [0, 1, 0])) ./ res_y;
    
    % update Dx, Dy, Dz
    Dx = mDx1 .* Dx + mDx2 .* CurlHx;
    Dy = mDy1 .* Dy + mDy2 .* CurlHy;
    Dz = Dz + mDz2 .* CurlHz + mDz3 .* (ICurlHz + 0.5 * CurlHz);
    
    % update Ex, Ey, Ez
    Ex = Dx ./ Eps_xx;
    Ey = Dy ./ Eps_yy;
    Ez = Dz ./ Eps_zz;
    
    % simple soft source (not recommended)
    % Ex(:, :, Ns_z) = Ex(:, :, Ns_z) + Ex_source(T);
    % Ey(:, :, Ns_z) = Ey(:, :, Ns_z) + Ey_source(T);
    
    % visualize
    if mod(T, 5) == 0
        clf;
        
        shift = 8;
        
        subplot(2, 3, 1);
        slice(Y, X, -Z, field_log_normalize(Ex, shift), 5, 5, 0);
        caxis([-shift, shift]);
        axis equal tight off;
        shading interp;
        hold on;
        
        plot3([5 5], [-5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z], '-k');
        hold on;
        plot3([5 5], [-5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z], '-k');
        hold on;
        surf(devX,devY,devZ,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8])
        hold on;
        plot3(dev_x,dev_y,dev_z,'k','LineWidth',1);
        hold on;
        plot3(dev_x',dev_y',dev_z','k','LineWidth',1);
        title('Ex');
         
        subplot(2, 3, 2);
        slice(Y, X, -Z, field_log_normalize(Ey, shift), 5, 5, 0);
        caxis([-shift, shift]);
        axis equal tight off;
        shading interp;
        hold on;
        plot3([5 5], [-5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z], '-k');
        hold on;
        plot3([5 5], [-5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z], '-k');
        hold on;        
        surf(devX,devY,devZ,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8])
        hold on;
        plot3(dev_x,dev_y,dev_z,'k','LineWidth',1);
        hold on;
        plot3(dev_x',dev_y',dev_z','k','LineWidth',1);        
        title('Ey');
        
        subplot(2, 3, 3);
        slice(Y, X, -Z, field_log_normalize(Ez, shift), 5, 5, 0);
        caxis([-shift, shift]);
        axis equal tight off;
        shading interp;
        hold on;
        plot3([5 5], [-5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z], '-k');
        hold on;
        plot3([5 5], [-5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z], '-k');
        hold on;        
        surf(devX,devY,devZ,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8])
        hold on;
        plot3(dev_x,dev_y,dev_z,'k','LineWidth',1);
        hold on;
        plot3(dev_x',dev_y',dev_z','k','LineWidth',1);
        title('Ez');
        
        subplot(2, 3, 4);
        slice(Y, X, -Z, field_log_normalize(Hx, shift), 5, 5, 0);
        caxis([-shift, shift]);
        axis equal tight off;
        shading interp;
        hold on;
        plot3([5 5], [-5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z], '-k');
        hold on;
        plot3([5 5], [-5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z], '-k');
        hold on;        
        surf(devX,devY,devZ,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8])
        hold on;
        plot3(dev_x,dev_y,dev_z,'k','LineWidth',1);
        hold on;
        plot3(dev_x',dev_y',dev_z','k','LineWidth',1);
        title('Hx');
        
        subplot(2, 3, 5);
        slice(Y, X, -Z, field_log_normalize(Hy, shift), 5, 5, 0);
        caxis([-shift, shift]);
        axis equal tight off;
        shading interp;
        hold on;
        plot3([5 5], [-5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z], '-k');
        hold on;
        plot3([5 5], [-5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z], '-k');
        hold on;        
        surf(devX,devY,devZ,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8])
        hold on;
        plot3(dev_x,dev_y,dev_z,'k','LineWidth',1);
        hold on;
        plot3(dev_x',dev_y',dev_z','k','LineWidth',1);        
        title('Hy');
        
        subplot(2, 3, 6);
        slice(Y, X, -Z, field_log_normalize(Hz, shift), 5, 5, 0);
        caxis([-shift, shift]);
        axis equal tight off;
        shading interp;
        hold on;
        plot3([5 5], [-5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [20.0 - N_layers_z1 * res_z 20.0 - N_layers_z1 * res_z], '-k');
        hold on;
        plot3([5 5], [-5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z],...
            '-k');
        hold on;
        plot3([-5 5], [5, 5],...
            [-20 + N_layers_z0 * res_z -20.0 + N_layers_z0 * res_z], '-k');
        hold on;        
        surf(devX,devY,devZ,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8])
        hold on;
        plot3(dev_x,dev_y,dev_z,'k','LineWidth',1);
        hold on;
        plot3(dev_x',dev_y',dev_z','k','LineWidth',1);
        title('Hz');
        
        t = Dt .*T;
        title_str = sprintf('3D FDTD Example, t = %.3e s', t);
        subtitle(title_str);
        
        F = getframe(gcf);
        writeVideo(vidObj, F);
        % pause(0.001);
    end
end
close(vidObj);