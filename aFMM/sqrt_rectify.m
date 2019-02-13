function LAM = sqrt_rectify(LAM)

LAM = sqrt(LAM);

[~, n] = size(LAM);

for i = 1 : n
    lam = LAM(i, i);
    lamreal = abs(real(lam));
    lamimag = abs(imag(lam));
    LAM(i, i) = lamreal + 1i * lamimag;
end

end