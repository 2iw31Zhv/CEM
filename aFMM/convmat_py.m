function PY = convmat_py(py, M, N, Lq)

py = fftshift(fft(py)) / Lq;

PY = zeros(M*N);
n0 = 1 + N;
n = (-floor(N/2):floor(N/2));

for nrow = 1 : N
    for mrow = 1 : M
        row = (nrow - 1) * M + mrow;
        for ncol = 1 : N
            for mcol = 1 : M
                col = (ncol - 1) * M + mcol;
                nfft = n(nrow) - n(ncol);
                if (mrow == mcol)
                    PY(row, col) = py(n0 + nfft);
                end
            end
        end
    end
end


end