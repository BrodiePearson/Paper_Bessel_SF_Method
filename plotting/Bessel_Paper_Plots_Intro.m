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
filedir_qg = "../simulations/output/data/final_simulation/";
filedir_copied_data = "../simulations/output/copied_data/SurfaceQG_n_512_visc_1.0e-17_order_8_hypovisc_0.005_order_-2_kf_7.0_F_1.0e-5_endtime_30000.0_beta_0.5/";

StrongBeta_2D = load("/Users/brodiepearson/GitHub/Paper_Anisotropic_2D_Structure_Functions/simulations/Output/Data/Equilibrated/Anisotropic2D_n_2048_drag_0.04_order_-2_visc_1.0e-21_order_8_kf_100.0_F_1.0e-5_betay_10.0_betax_0.0_3000.mat");
SQG = load(filedir_copied_data+"SurfaceQG_n_512_visc_1.0e-17_order_8_hypovisc_0.005_order_-2_kf_7.0_F_1.0e-5_endtime_30000.0_beta_0.5_30040.mat");
QG_mid = load(filedir_qg+"multilayerqg_2layer_512_beta_5.0_Ld_0.35_75.mat");
QG = load(filedir_qg+"multilayerqg_2layer_512_beta_5.0_Ld_0.35_100.mat");

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

%% Create a QG plot showing evolution of energy through time

experiment_name = "multilayerqg_2layer_512_beta_5.0_Ld_0.35";

filedir = "../simulations/output/data/final_simulation/";

files = dir(filedir + experiment_name +"_*.mat");

% Sort the snapshots from earliest to latest
temp1 = struct2table(files); % convert the struct array to a table
temp2 = sortrows(temp1, 'datenum'); % sort the table by 'DOB'
files = table2struct(temp2);
file_no = length(files);

for layer=1:2
    for nn=1:file_no
        S = load(filedir + files(nn).name);
        kinetic_energy(layer,nn) = S.KE(layer);
        time(layer,nn)=files(nn).datenum;
        potential_energy(nn) = S.PE;
        vert_flux(nn) = S.Vert_flux;
    end
end

h_qgmaps=figure(21)
set(h_qgmaps,'Position',[10 10 1600 400])
ax(1) = subplot(1,3,1)
upper_layer = plot(kinetic_energy(1,:)/0.2, 'k','LineWidth',2);
hold on
lower_layer = plot(kinetic_energy(2,:)/0.8, 'r','LineWidth',2);
Vert_Flux = plot(vert_flux(:), 'r','LineWidth',2);
%PE = plot(potential_energy(:)*10, 'b:','LineWidth',2);
plot ([20 20], [0 3000],'k:','LineWidth',2)
plot ([75 75], [0 3000],'k--','LineWidth',2)
plot ([100 100], [0 3000],'k--','LineWidth',2)
ylabel('Energy per unit depth (m s^{-2})')
xlabel('Time (s)')
legend([upper_layer lower_layer], ...
    'Upper layer KE','Lower Layer KE', ...
    'Location','North');
xlim([0 105])
set(gca,'fontsize', 14);

ax(2) = subplot(1,3,2)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG_mid.q(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('Time = 75')
set(gca,'fontsize', 14);
colormap(ax(2), color_bar_qgpv)
hcb=colorbar
caxis([-3e3 3e3])
hcb.Label.String = 'QG Potential Vorticity';
hcb.Location = 'southoutside';

ax(3) = subplot(1,3,3)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG.q(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('Time = 100')
set(gca,'fontsize', 14);
colormap(ax(3), color_bar_qgpv)
hcb=colorbar
caxis([-3e3 3e3])
hcb.Label.String = 'QG Potential Vorticity';
hcb.Location = 'southoutside';

%% Alt panel

QG_midalt = load(filedir_qg+"multilayerqg_2layer_512_beta_5.0_Ld_0.35_41.mat");

h_qgmaps=figure(30)
set(h_qgmaps,'Position',[10 10 1600 400])
ax(1) = subplot(1,3,1)
upper_layer = plot(kinetic_energy(1,:)/0.2, 'k','LineWidth',2);
hold on
lower_layer = plot(kinetic_energy(2,:)/0.8, 'r','LineWidth',2);
%Vert_Flux = plot(vert_flux(:), 'r','LineWidth',2);
%PE = plot(potential_energy(:)*10, 'b:','LineWidth',2);
plot ([20 20], [0 3000],'k:','LineWidth',2)
plot ([41 41], [0 3000],'k--','LineWidth',2)
plot ([100 100], [0 3000],'k--','LineWidth',2)
ylabel('Energy per unit depth (m s^{-2})')
xlabel('Time (s)')
legend([upper_layer lower_layer], ...
    'Upper layer KE','Lower Layer KE', ...
    'Location','North');
xlim([0 105])
set(gca,'fontsize', 14);

ax(2) = subplot(1,3,2)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG_midalt.q(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('Time = 75')
set(gca,'fontsize', 14);
colormap(ax(2), color_bar_qgpv)
hcb=colorbar
caxis([-3e3 3e3])
hcb.Label.String = 'QG Potential Vorticity';
hcb.Location = 'southoutside';

ax(3) = subplot(1,3,3)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG.q(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('Time = 100')
set(gca,'fontsize', 14);
colormap(ax(3), color_bar_qgpv)
hcb=colorbar
caxis([-3e3 3e3])
hcb.Label.String = 'QG Potential Vorticity';
hcb.Location = 'southoutside';





%% Create a QG plot showing evolution of energy through time

h_qgmaps=figure(35)
set(h_qgmaps,'Position',[10 10 1600 400])


QG_temp1 = load(filedir_qg+"multilayerqg_2layer_512_beta_5.0_Ld_0.35_41.mat");
QG_temp2 = load(filedir_qg+"multilayerqg_2layer_512_beta_5.0_Ld_0.35_42.mat");
QG_temp3 = load(filedir_qg+"multilayerqg_2layer_512_beta_5.0_Ld_0.35_43.mat");

ax(1) = subplot(1,3,1)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG_temp1.q(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('Time = 100')
set(gca,'fontsize', 14);
colormap(ax(1), color_bar_qgpv)
hcb=colorbar
caxis([-3e3 3e3])
hcb.Label.String = 'QG Potential Vorticity';
hcb.Location = 'southoutside';

ax(2) = subplot(1,3,2)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG_temp2.q(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('Time = 75')
set(gca,'fontsize', 14);
colormap(ax(2), color_bar_qgpv)
hcb=colorbar
caxis([-3e3 3e3])
hcb.Label.String = 'QG Potential Vorticity';
hcb.Location = 'southoutside';

ax(3) = subplot(1,3,3)
imagesc(x_grid_QG(1:end),y_grid_QG(1:end), QG_temp3.q(1:end,1:end,1)')
ylim([0 2*pi])
xlim([0 2*pi])
xticks([0, pi, 2*pi])
xticklabels({'0' '\pi' '2\pi'})
yticks([0, pi, 2*pi])
yticklabels({'0', '\pi', '2\pi'})
%xlabel('x (m)')
%ylabel('y (m)')
title('Time = 100')
set(gca,'fontsize', 14);
colormap(ax(3), color_bar_qgpv)
hcb=colorbar
caxis([-3e3 3e3])
hcb.Label.String = 'QG Potential Vorticity';
hcb.Location = 'southoutside';



