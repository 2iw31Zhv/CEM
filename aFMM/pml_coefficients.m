function p = pml_coefficients(c_PML, Lq, Dq, M)

kappa = atan(pi/2);
q = linspace(-0.5 * Lq, 0.5 * Lq, M);
d_beta_q = zeros(M, 1);

for i = 1 : M
    if abs(q(i)) <= 0.5 * Dq
        d_beta_q(i) = 0.0;
    elseif q(i) > 0.5 * Dq
        d_beta_q(i) = - 2 * kappa * c_PML / (Lq - Dq)...
            * (tan(tan(kappa * (2 * abs(q(i)) - Dq) / (Lq - Dq))))^2 ...
         / (cos(kappa * (2 * abs(q(i)) - Dq) / (Lq - Dq)))^2;
    else
        d_beta_q(i) = - 2 * kappa * c_PML / (Lq - Dq)...
            * (tan(tan(kappa * (2 * abs(q(i)) - Dq) / (Lq - Dq))))^2 ...
         / (cos(kappa * (2 * abs(q(i)) - Dq) / (Lq - Dq)))^2;
    end
end

inv_d_T = 1 ./ (1 + 1i .* d_beta_q);
p = fftshift(fftn(inv_d_T)) / Lq;

end