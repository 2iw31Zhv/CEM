close all;
clear variables;

load s4matrices.mat
load k.mat
load rcwa_sg.mat

Srect = sk2srect(S13, k);

% fig1 = figure('Color', 'w');
% imagesc(abs(Srect.S11(1:121, 1:121)));
% caxis([0, 2.5]);
% colorbar;
% 
% fig2 = figure('Color', 'w');
% imagesc(abs(SG.S11(1:121, 1:121)));
% caxis([0, 2.5]);
% colorbar;
% 
% saveas(fig1, 'Srect_S11_xx.png');
% saveas(fig2, 'SG_S11_xx.png');

fig3 = figure('Color', 'w');
imagesc(abs(Srect.S12(122:242, 122:242)));
caxis([0, 2.5]);
colorbar;

fig4 = figure('Color', 'w');
imagesc(abs(SG.S12(122:242, 122:242)));
caxis([0, 2.5]);
colorbar;

saveas(fig3, 'Srect_S12_yy.png');
saveas(fig4, 'SG_S12_yy.png');
