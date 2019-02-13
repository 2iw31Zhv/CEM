function [ER, UR] = set_circle_pattern(Nx, Ny, Cx, Cy, Lx, Ly, r,...
    eps_r, eps0, mu_r, mu0)
dx = Lx / Nx;
dy = Ly / Ny;
ER = ones(Nx, Ny);
UR = ones(Nx, Ny);

for i  = 1 : Nx
    for j = 1 : Ny
        x = Cx - 0.5 * Lx + (i - 0.5) * dx;
        y = Cy - 0.5 * Ly + (j - 0.5) * dy;
        
        if ((x - Cx)^2 + (y - Cy)^2 <= r^2)
            ER(i, j) = eps_r;
            UR(i, j) = mu_r;
        else
            ER(i, j) = eps0;
            UR(i, j) = mu0;
        end
    end
end

end