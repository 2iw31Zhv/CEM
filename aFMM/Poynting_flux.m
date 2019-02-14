function POWER = Poynting_flux(c, W, V, eta0)
n = size(c);
n = n / 2;

eT = W * c;
ex = eT(1:n);
ey = eT((n+1):2*n);

hT = - eta0 * V * c;
hx = hT(1:n);
hy = hT((n+1):2*n);

POWER = -0.5 * real(sum(ex .* conj(hy) - ey .* conj(hx)));

end