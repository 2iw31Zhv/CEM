% =========================================================================
% 2D Plane Wave Expansion Method Example (E Mode)
% =========================================================================
% Author: Ziwei Zhu
% Reference: http://emlab.utep.edu/ee5390cem.htm
% =========================================================================

close all;
clear variables;

% solve region
% define a x a square region with periodic condition
a = 1.0; % m
% samples of the Fourier basis
% remember we assume them to be odd
P = 7;
Q = 7;

% devive parameter
% a round hole with radius r
r = 0.35 * a;
% relative dielectric constant of substrate
epsilon_r_out = 1.5;
epsilon_r_in = 2.5;

% physical constant

% set grid parameters
Nx = 512;
Ny = 512;
Dx = a / Nx;
Dy = a / Ny;


% define the permittivity and permeability field
UR = ones(Nx, Ny);
ER = epsilon_r_in * ones(Nx, Ny);

for i = 1 : Nx
    for j = 1 : Ny
        real_x = -0.5 * a + Dx * (i - 0.5);
        real_y = -0.5 * a + Dy * (j - 0.5);
        
        if (real_x * real_x + real_y * real_y > r * r)
            ER(i, j) = epsilon_r_out;
        end
    end
end

% convert permittivity and permeability to convolution matrices
ERC = convmat(ER, P, Q);
URC = convmat(UR, P, Q);

%imagesc(1:P, 1:Q, abs(ERC));
%axis equal tight;

% compute the key points in the Irreducible Brillouin Zone
T1 = 2 * pi / a * [1; 0];
T2 = 2 * pi / a * [0; 1];
Gamma = [0; 0];
CHI = 0.5 * T1;
M = 0.5 * T1 + 0.5 * T2;

Nb = 200;
NB = 1 + 2  * Nb + round(sqrt(2) * Nb);
beta = zeros(2, NB);

% interpolation
for i = 1 : NB
    if i <= (Nb + 1)
        alpha = (i-1) / Nb;
        beta(:, i) = (1 - alpha) * Gamma + alpha * CHI;
    elseif i <= (2*Nb+1)
        alpha = (i - Nb - 1) / Nb;
        beta(:, i) = (1 - alpha) * CHI + alpha * M;
    else
        alpha = (i - 2*Nb - 1) / round(sqrt(2) * Nb);
        beta(:, i) = (1 - alpha) * M + alpha * Gamma;
    end
end

energy_band = zeros(P*Q, NB);

for i = 1 : NB
    
    bx = beta(1, i);
    by = beta(2, i);
    
    p = (-floor(P/2):floor(P/2));
    q = (-floor(Q/2):floor(Q/2));
    
    kx = bx - 2*pi/a * p;
    ky = by - 2*pi/a * q;
    
    [ky, kx] = meshgrid(ky, kx);
    KX = diag(sparse(kx(:)));
    KY = diag(sparse(ky(:)));
    
    A = KX * URC^(-1) * KX + KY * URC^(-1) * KY;
    B = ERC;
    
    % solve the general eigenvalue problem
    [V, D] = eig(A, B);
    energy_band(:, i) = sqrt(diag(D)) * a / 2 / pi;
end


fig = figure('Color', 'w', 'Position', [100 300 800 400]);

for i = 1 : (P * Q)
    scatter(1:NB, energy_band(i, :), 1);
    hold on;
end

axis([1, NB, 0, 1]);
ylabel('normalized frequency \omega a / 2 \pi c');
set(gca,'xtick',[1 Nb+1 2*Nb+1 NB]);
set(gca,'xticklabel',{'\Gamma','X','M','\Gamma'});
xlabel('Bloch wave vector \beta');
title('Energy band diagram of round waveguide');
box on;
saveas(fig, 'energyband.png');

% the final recorded would be Gamma
% the first with non-zero energy
X = -0.5 * a + Dx * ((1:Nx) - 0.5);
Y = -0.5 * a + Dy * ((1:Ny) - 0.5);
[Y, X] = meshgrid(Y, X);

for i = 1 : 8
    mode_visualize(V, i, [0;0], P, Q, X, Y, 0.35);
end

