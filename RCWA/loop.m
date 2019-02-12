close all;
clear variables;

lambda = linspace(2.0, 5.0, 100);
[~, n] = size(lambda);

REF = zeros(n, 1);
TRN = zeros(n, 1);
SUM = zeros(n, 1);

for i = 1 : n
    lambda0 = lambda(i);
    [REF(i), TRN(i), SUM(i)] = rcwa3d(lambda0);
end

fig = figure('Color', 'w');
plot(lambda, REF, '-r', 'LineWidth', 3);
hold on;
plot(lambda, TRN, '-g', 'LineWidth', 2);
hold on;
plot(lambda, SUM, '-k', 'LineWidth', 1);
axis([2.0, 5.0, 0, 1.2]);
legend('REF', 'TRN', 'Conservation');
saveas(fig, 'images/conservation.png');

