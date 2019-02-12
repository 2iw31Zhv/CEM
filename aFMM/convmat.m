function C = convmat(A, P, Q, R)
% Explanation:
% convert a permittivity/permeability matrix from real space to a
% convolutional matrix in Fourier space

% Input:
% C = convmat(A, P)                 for 1D problem
% C = convmat(A, P, Q)              for 2D problem
% C = convmat(A, P, Q, R)           for 3D problem


[Nx, Ny, Nz] = size(A);

if nargin == 2
    Q = 1;
    R = 1;
elseif nargin == 3
    R = 1;
end

NH = P * Q * R;
p = (-floor(P/2):floor(P/2));
q = (-floor(Q/2):floor(Q/2));
r = (-floor(R/2):floor(R/2));

A = fftshift(fftn(A)) / (Nx * Ny * Nz);

p0 = 1 + floor(Nx/2);
q0 = 1 + floor(Ny/2);
r0 = 1 + floor(Nz/2);

C = zeros(NH, NH);

for rrow = 1 : R
    for qrow = 1 : Q
        for prow = 1 : P
            row = (rrow - 1) * P * Q + (qrow - 1) * P + prow;
            for rcol = 1 : R
                for qcol = 1 : Q
                    for pcol = 1 : P
                        col = (rcol - 1) * P * Q + (qcol-1) * P + pcol;
                        pfft = p(prow) - p(pcol);
                        qfft = q(qrow) - q(qcol);
                        rfft = r(rrow) - r(rcol);
                        C(row, col) = A(p0 + pfft, q0 + qfft, r0 + rfft);
                    end
                end
            end
        end
    end
end

end

