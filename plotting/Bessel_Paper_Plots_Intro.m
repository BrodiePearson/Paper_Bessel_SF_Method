%% Bessel function tests

clear all

size_of_font = 20;

%% Plot spectral fluxes

figure(99)
z = 0:0.1:20;
J = zeros(5,201);
xJ = zeros(5,201);
testJ = zeros(1,201);
for i = 0:3
    J(i+1,:) = besselj(i,z);
    zJ(i+1,:) = z.*besselj(i,z);
end
plot(z,J(2,:), 'k')
hold on
plot(z,J(3,:), 'k--')
plot(z,J(4,:), 'k:')
grid on
legend('J_1','J_2','J_3','Location','Best')
xlabel('z')
ylabel('J_{\nu}(z)')
set(gca,'fontsize', size_of_font);

%% Plot vorticity fields from snapshots of each simulation

addpath ../../AR6_Chapter9/FGD/Matlab_Functions/

color_bar_vorticity = IPCC_Get_Colorbar('temperature_d', 21, false);
color_bar_buoyancy = IPCC_Get_Colorbar('precip_d', 21, false);
color_bar_qgpv = IPCC_Get_Colorbar('wind_d', 21, false);

filedir = "../simulations/output/data/";
filedir_copied_data = "../simulations/output/copied_data/SurfaceQG_n_512_visc_1.0e-17_order_8_hypovisc_0.005_order_-2_kf_7.0_F_1.0e-5_endtime_30000.0_beta_0.5/";

StrongBeta_2D = load("/Users/brodiepearson/GitHub/Paper_Anisotropic_2D_Structure_Functions/simulations/Output/Data/Equilibrated/Anisotropic2D_n_2048_drag_0.04_order_-2_visc_1.0e-21_order_8_kf_100.0_F_1.0e-5_betay_10.0_betax_0.0_3000.mat");
SQG = load(filedir_copied_data+"SurfaceQG_n_512_visc_1.0e-17_order_8_hypovisc_0.005_order_-2_kf_7.0_F_1.0e-5_endtime_30000.0_beta_0.5_30040.mat");
QG = load(filedir+"multilayerqg_2layer_512_beta_5.0_Ld_0.35_46.mat");

x_grid_2D = StrongBeta_2D.dx:StrongBeta_2D.dx:StrongBeta_2D.Lx;
y_grid_2D = StrongBeta_2D.dy:StrongBeta_2D.dy:StrongBeta_2D.Ly;
x_grid_SQG = SQG.dx:SQG.dx:SQG.Lx;
y_grid_SQG = SQG.dy:SQG.dy:SQG.Ly;
x_grid_QG = QG.dx:QG.dx:QG.Lx;
y_grid_QG = QG.dy:QG.dy:QG.Ly;

hmaps=figure(20)
set(hmaps,'Position',[10 10 1600 600])
ax(1) = subplot(2,3,1)
imagesc(x_grid_2D(1:end),y_grid_2D(1:end), StrongBeta_2D.zeta(1:end,1:end)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
xlabel('x (m)')
ylabel('y (m)')
title('2D Turbulence')
set(gca,'fontsize', 20);
colormap(ax(1), color_bar_vorticity)
hcb=colorbar
caxis([-5 5])
hcb.Label.String = 'Vorticity [s^{-1}]';
hcb.Location = 'southoutside';

ax(2) = subplot(2,3,2)
imagesc(x_grid_SQG(1:end),y_grid_SQG(1:end), SQG.b(1:end,1:end)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title(' Surface Quasigeostrophic (SQG) Turbulence')
set(gca,'fontsize', 20);
colormap(ax(2), color_bar_buoyancy)
hcb=colorbar
caxis([-0.1 0.1])
hcb.Label.String = 'Buoyancy [s^{-1}]';
hcb.Location = 'southoutside';

ax(3) = subplot(2,3,3)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG.q(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('Quasi-geostrophic (QG) Turbulence')
set(gca,'fontsize', 20);
colormap(ax(3), color_bar_qgpv)
hcb=colorbar
caxis([-1e3 1e3])
hcb.Label.String = 'QG Potential Vorticity [s^{-1}]';
hcb.Location = 'southoutside';