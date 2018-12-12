function S = redheffer_star_product(SA, SB)
% this function combines two scatter matrices into one total scatter matrix
% input: S1, S2 two partial scatter matrices
% output: S the total scatter matrix

[~, n] = size(SA); % assume SA and SB are of the same size
n = n / 2;

SA11 = SA(1:n, 1:n);
SA12 = SA(1:n, (n+1):2*n);
SA21 = SA((n+1):2*n, 1:n);
SA22 = SA((n+1):2*n, (n+1):2*n);

SB11 = SB(1:n, 1:n);
SB12 = SB(1:n, (n+1):2*n);
SB21 = SB((n+1):2*n, 1:n);
SB22 = SB((n+1):2*n, (n+1):2*n);

S = zeros(2*n, 2*n);
I = eye(n);

S(1:n, 1:n) = SA11 + SA12 * (I - SB11 * SA22)^(-1) * SB11 * SA21;
S(1:n, (n+1):2*n) = SA12 * (I - SB11 * SA22)^(-1) * SB12;
S((n+1):2*n, 1:n) = SB21 * (I - SA22 * SB11)^(-1) * SA21;
S((n+1):2*n, (n+1):2*n) = SB22 + SB21 * (I - SA22 * SB11)^(-1) * SA22 * SB12;

end