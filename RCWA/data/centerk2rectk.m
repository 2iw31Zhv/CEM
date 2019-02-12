function Srect = centerk2rectk(S, K, M, N)

% =========================================================================
% input : S, a scatter matrix with k vectors specified by K
% output: Srect, a scatter matrix, with M * N k vectors sampled
% homogenously in a rectangle
% =========================================================================

Srect = zeros(M*N, M*N);

[nS, ~] = size(S);

for i = 1 : nS
    for j = 1 : nS
        kxi = K(i, 1);
        kyi = K(i, 2);
        kxj = K(j, 1);
        kyj = K(j, 2);
        
        mi = round(kxi / 2 / pi) + floor(M/2) + 1;
        ni = round(kyi / 2 / pi) + floor(N/2) + 1;   
        mj = round(kxj / 2 / pi) + floor(M/2) + 1;
        nj = round(kyj / 2 / pi) + floor(N/2) + 1;
        
        
        Srect(mi + (ni-1) * M, mj + (nj-1) * M) = S(i, j);
    end
end


end