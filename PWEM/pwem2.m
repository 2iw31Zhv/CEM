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
P = 7;
Q = 7;

% devive parameter
% a round hole with radius r
r = 0.35 * a;
% relative dielectric constant of substrate
epsilon_r = 9.0;

% physical constant

% set grid parameters
Nx = 512;
Ny = 512;
Dx = a / Nx;
Dy = a / Ny;


% define the permittivity and permeability field
UR = ones(Nx, Ny);
ER = ones(Nx, Ny);

for i = 1 : Nx
    for j = 1 : Ny
        real_x = -0.5 * a + Dx * (i - 0.5);
        real_y = -0.5 * a + Dy * (j - 0.5);
        
        if (real_x * real_x + real_y * real_y > r * r)
            ER(i, j) = epsilon_r;
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

Nb = 20;
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

for i = 1 : NB
    
end