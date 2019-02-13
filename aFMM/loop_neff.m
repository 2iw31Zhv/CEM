close all;
clear variables;

n = 500;
freq_min = 0.01;
freq_max = 5.01;
freq = linspace(freq_min, freq_max, n);

neff = zeros(n, 1);

for i = 1 : n
    freq0 = freq(i);
    neff(i) = aFMM(freq0);
end

fig = figure('Color', 'w');
plot(freq, neff, '-k', 'LineWidth', 1);
axis([freq_min, freq_max, 0, 5.0]);
saveas(fig, 'neff_radius_0.4.png');
