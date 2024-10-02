%% Bessel function paper plots (maps of fields from all simulations)

clear all

%% Plot vorticity fields from three beta-effect simulations

addpath ../../AR6_Chapter9/FGD/Matlab_Functions/

color_bar_vort = IPCC_Get_Colorbar('temperature_d', 21, false);
color_bar_b = IPCC_Get_Colorbar('precip_d', 21, false);

% SB = load("/Users/brodiepearson/GitHub/Paper_Anisotropic_2D_Structure_Functions/simulations/Output/Data/Equilibrated/Anisotropic2D_n_2048_drag_0.04_order_-2_visc_1.0e-21_order_8_kf_100.0_F_1.0e-5_betay_10.0_betax_0.0_3000.mat");
% SBd = load("/Users/brodiepearson/GitHub/Paper_Anisotropic_2D_Structure_Functions/simulations/Output/Data/Equilibrated/Anisotropic2D_n_2048_drag_0.04_order_-2_visc_1.0e-21_order_8_kf_100.0_F_1.0e-5_betay_7.071067811865475_betax_7.071067811865475_3000.mat");
% WB = load("/Users/brodiepearson/GitHub/Paper_Anisotropic_2D_Structure_Functions/simulations/Output/Data/Equilibrated/Anisotropic2D_n_2048_drag_0.04_order_-2_visc_1.0e-21_order_8_kf_100.0_F_1.0e-5_betay_0.7071067811865475_betax_0.7071067811865475_3000.mat");
% SBn = load("/Users/brodiepearson/GitHub/Paper_Anisotropic_2D_Structure_Functions/simulations/Output/Data/Equilibrated/Anisotropic2D_n_2048_drag_0.04_order_-2_visc_1.0e-21_order_8_kf_40.0_F_1.0e-5_betay_10.0_betax_0.0_3000.mat");

data_path = "../simulations/output/copied_data/"
TwoD = load(data_path+"2D_Equilibrated/Anisotropic2D_n_2048_drag_0.04_order_-2_visc_1.0e-21_order_8_kf_100.0_F_1.0e-5_betay_10.0_betax_0.0_3000.mat");
SQG = load(data_path+"SurfaceQG_n_512_visc_1.0e-17_order_8_hypovisc_0.005_order_-2_kf_7.0_F_1.0e-5_endtime_30000.0_beta_0.5/SurfaceQG_n_512_visc_1.0e-17_order_8_hypovisc_0.005_order_-2_kf_7.0_F_1.0e-5_endtime_30000.0_beta_0.5_30040.mat");
QG = load("/Users/brodiepearson/GitHub/Paper_Bessel_SF_Method/simulations/output/data/final_simulation/multilayerqg_2layer_512_beta_5.0_Ld_0.35_100.mat");
%SBn = load("/Users/brodiepearson/GitHub/Paper_Anisotropic_2D_Structure_Functions/simulations/Output/Data/Equilibrated/Anisotropic2D_n_2048_drag_0.04_order_-2_visc_1.0e-21_order_8_kf_40.0_F_1.0e-5_betay_10.0_betax_0.0_3000.mat");


x_grid_2D = TwoD.dx:TwoD.dx:TwoD.Lx;
y_grid_2D = TwoD.dy:TwoD.dy:TwoD.Ly;

x_grid_SQG = SQG.dx:SQG.dx:SQG.Lx;
y_grid_SQG = SQG.dy:SQG.dy:SQG.Ly;

x_grid_QG = QG.dx:QG.dx:QG.Lx;
y_grid_QG = QG.dy:QG.dy:QG.Ly;

hmaps=figure(20)
set(hmaps,'Position',[10 10 1300 1300])

ax(1) = subplot(3,2,1)
imagesc(x_grid_2D(1:end),y_grid_2D(1:end), TwoD.zeta(1:end,1:end)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
ylabel('y (m)')
title('2D Turbulence (\omega)')
set(gca,'fontsize', 20);
caxis([-5 5])
colormap(ax(1), color_bar_vort)
hcb_vort=colorbar
hcb_vort.Label.String = 'Vorticity \omega [s^{-1}], q_1 [10 s^{-1}], and q_2 [s^{-1}]';
hcb_vort.Location = 'southoutside';


ax(3) = subplot(3,2,3)
imagesc(x_grid_SQG(1:end),y_grid_SQG(1:end), SQG.b(1:end,1:end)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
xlabel('x (m)')
%ylabel('y (m)')
title(' SQG Turbulence (b)')
set(gca,'fontsize', 20);
colormap(ax(3), color_bar_b)
caxis([-0.25 0.25])
hcb_b=colorbar
hcb_b.Label.String = 'Buoyancy [m^2 s^{-1}]';
hcb_b.Location = 'southoutside';


ax(2) = subplot(3,2,2)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG.q(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('QG Upper Layer (q_1)')
set(gca,'fontsize', 20);
colormap(ax(2), color_bar_vort)
caxis([-50 50])

ax(4) = subplot(3,2,4)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG.q(1:end,1:end,2)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
xlabel('x (m)')
%ylabel('y (m)')
title('QG Lower Layer (q_2)')
set(gca,'fontsize', 20);
colormap(ax(4), color_bar_vort)
caxis([-5 5])

%% QG Streamfunction map

figure(11)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG.psi(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])

yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('QG Upper Layer Streamfunction (\psi_1)')
set(gca,'fontsize', 20);
colormap(color_bar_vort)
colorbar
%caxis([-50 50])


%% QG Streamfunction map

figure(12)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG.psi(1:end,1:end,2)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])

yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('QG Lower Layer Streamfunction (\psi_1)')
set(gca,'fontsize', 20);
colormap(color_bar_vort)
colorbar
%caxis([-50 50])

%% QG velocity map

figure(13)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG.u(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])

yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('QG Upper Layer zonal velocity (u_1)')
set(gca,'fontsize', 20);
colormap(color_bar_vort)
colorbar
%caxis([-50 50])

%% QG velocity map

figure(14)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG.v(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])

yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('QG Upper Layer meridional velocity (u_1)')
set(gca,'fontsize', 20);
colormap(color_bar_vort)
colorbar