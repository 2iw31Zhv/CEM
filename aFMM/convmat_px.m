function PX = convmat_px(px, M, N)

PX = zeros(M*N);
m0 = 1 + M;
m = (-floor(M/2):floor(M/2));

for nrow = 1 : N
    for mrow = 1 : M
        row = (nrow - 1) * M + mrow;
        for ncol = 1 : N
            for mcol = 1 : M
                col = (ncol - 1) * M + mcol;
                mfft = m(mrow) - m(mcol);
                if (nrow == ncol)
                    PX(row, col) = px(m0 + mfft);
                end
            end
        end
    end
end

end