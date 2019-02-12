function Srect = sk2srect(S, K)

absK = sqrt((K(:, 1) ./ 2 ./ pi).^2 + (K(:, 2) ./ 2 ./ pi).^2);
M = 2 * round(max(absK)) + 1;
N = M;

n = M * N;
n2 = 2 * n;

Srect = struct;
Srect.S11 = zeros(n2, n2);
Srect.S12 = zeros(n2, n2);
Srect.S21 = zeros(n2, n2);
Srect.S22 = zeros(n2, n2);

[nS, ~] = size(S);
nS2 = nS / 2;
nS = nS / 4;

S11 = S(1:nS2, 1:nS2);
S12 = S(1:nS2, (nS2+1):2*nS2);
S21 = S((nS2+1):2*nS2, 1:nS2);
S22 = S((nS2+1):2*nS2, (nS2+1):2*nS2);

Srect.S11(1:n, 1:n) = centerk2rectk(S11(1:nS, 1:nS), K, M, N);
Srect.S11(1:n, (n+1):2*n)  = centerk2rectk(S11(1:nS, (nS+1):2*nS), K, M, N);
Srect.S11((n+1):2*n, 1:n)  = centerk2rectk(S11((nS+1):2*nS, 1:nS), K, M, N);
Srect.S11((n+1):2*n, (n+1):2*n)  = centerk2rectk(S11((nS+1):2*nS, (nS+1):2*nS), K, M, N);

Srect.S12(1:n, 1:n) = centerk2rectk(S12(1:nS, 1:nS), K, M, N);
Srect.S12(1:n, (n+1):2*n)  = centerk2rectk(S12(1:nS, (nS+1):2*nS), K, M, N);
Srect.S12((n+1):2*n, 1:n)  = centerk2rectk(S12((nS+1):2*nS, 1:nS), K, M, N);
Srect.S12((n+1):2*n, (n+1):2*n)  = centerk2rectk(S12((nS+1):2*nS, (nS+1):2*nS), K, M, N);

Srect.S21(1:n, 1:n) = centerk2rectk(S21(1:nS, 1:nS), K, M, N);
Srect.S21(1:n, (n+1):2*n)  = centerk2rectk(S21(1:nS, (nS+1):2*nS), K, M, N);
Srect.S21((n+1):2*n, 1:n)  = centerk2rectk(S21((nS+1):2*nS, 1:nS), K, M, N);
Srect.S21((n+1):2*n, (n+1):2*n)  = centerk2rectk(S21((nS+1):2*nS, (nS+1):2*nS), K, M, N);

Srect.S22(1:n, 1:n) = centerk2rectk(S22(1:nS, 1:nS), K, M, N);
Srect.S22(1:n, (n+1):2*n)  = centerk2rectk(S22(1:nS, (nS+1):2*nS), K, M, N);
Srect.S22((n+1):2*n, 1:n)  = centerk2rectk(S22((nS+1):2*nS, 1:nS), K, M, N);
Srect.S22((n+1):2*n, (n+1):2*n)  = centerk2rectk(S22((nS+1):2*nS, (nS+1):2*nS), K, M, N);

end