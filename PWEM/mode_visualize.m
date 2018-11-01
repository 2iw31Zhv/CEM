function mode_visualize(V, i, beta, P, Q, X, Y, r)
    close all;
    [Nx, Ny] = size(X);
    S = reshape(V(:, i), [P, Q]);
       
    nxc = ceil(Nx/2);
    nx1 = nxc - floor(P/2);
    nx2 = nxc + floor(P/2);
    nyc = ceil(Ny/2);
    ny1 = nyc - floor(Q/2);
    ny2 = nyc + floor(Q/2);
    
    
    sf = zeros(Nx,Ny);
    sf(nx1:nx2,ny1:ny2) = S;
    az = ifft2(ifftshift(sf));
    az = az / max(abs(az(:)));
    
    phase = exp(-1i*( beta(1) * X + beta(2) * Y));
    
    imagesc(abs(phase .* az));
    axis equal tight;
    
    hold on;
    phi = linspace(0,2*pi,50);
    x = Nx / 2 + r * Nx * cos(phi);
    y = Ny / 2 + r * Ny * sin(phi);
    plot(x, y, '-k');
    
    axis equal tight;
    title_str = sprintf('mode_%d.png', i);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    saveas(gcf, title_str);
end