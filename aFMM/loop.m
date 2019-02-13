close all;
clear variables;

freq = linspace(0.01, 5.01, 500);
[~, n] = size(freq);

REF = zeros(n, 1);
TRN = zeros(n, 1);
SUM = zeros(n, 1);

for i = 1 : n
    freq0 = freq(i);
    [REF(i), TRN(i), SUM(i)] = aFMM(freq0);
end

fig = figure('Color', 'w');
plot(freq, REF, '-r', 'LineWidth', 3);
hold on;
plot(freq, TRN, '-g', 'LineWidth', 2);
hold on;
plot(freq, SUM, '-k', 'LineWidth', 1);
axis([0.0, 5.0, 0, 3]);
legend('REF', 'TRN', 'Conservation');
saveas(fig, 'conservation.png');

