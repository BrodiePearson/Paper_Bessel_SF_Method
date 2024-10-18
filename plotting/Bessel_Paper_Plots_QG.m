%% Bessel function tests

clear all

size_of_font = 20;
legend_font_size = 14;
J_1_maximum = 2;
J_2_maximum = 3.1;
J_3_maximum = 4.3;

xmin = 0.7;
xmax = 250;

Rossby_deformation_radius=0.35;
Rossby_deformation_k = 2*pi/Rossby_deformation_radius;

experiment_name = "multilayerqg_2layer_512_beta_5.0_Ld_"+Rossby_deformation_radius;

quantiles = 0.75;

no_of_exps = size(experiment_name,2);

k_Enstrophy_Cascade_End = 4e2;
r_Enstrophy_Cascade_End = 2*pi/k_Enstrophy_Cascade_End;

k_forcing_estimate = 1e1;
r_forcing_estimate = 2*pi/k_forcing_estimate;
path = "../analysis/processed_data/";

Diagonal = load(path + 'Structure_Functions_'+experiment_name+'_Diagonal.mat');
Meridional = load(path + 'Structure_Functions_'+experiment_name+'_Meridional.mat');
Zonal = load(path + 'Structure_Functions_'+experiment_name+'_Zonal.mat');
OffDiagonal = load(path + 'Structure_Functions_'+experiment_name+'_Off-Diagonal.mat');


%% Plot spectral fluxes

Spectral_Flux = load(path + 'Spectral_Fluxes_'+experiment_name+'.mat');

Patch_k_Bounds = [Spectral_Flux.Wavenumber' fliplr(Spectral_Flux.Wavenumber')];
Patch_k_Bounds_Enstrophy = [k_forcing_estimate, max(Spectral_Flux.Wavenumber), ...
    max(Spectral_Flux.Wavenumber), k_forcing_estimate];
%Patch_k_Bounds_Energy = [k_Energy_Cascade_End, k_forcing_estimate, ...
%    k_forcing_estimate, k_Energy_Cascade_End];

Spectral_Flux.qq_Flux_max = (quantile(Spectral_Flux.qq_Flux_snapshots ...
    ,quantiles,3))-Spectral_Flux.qq_Flux(1);
Spectral_Flux.qq_Flux_min = (quantile(Spectral_Flux.qq_Flux_snapshots ...
    ,1-quantiles,3))-Spectral_Flux.qq_Flux(1);
Spectral_Flux.KE_Flux_max = (quantile(Spectral_Flux.KE_Flux_snapshots ...
    ,quantiles,3));
Spectral_Flux.KE_Flux_min = (quantile(Spectral_Flux.KE_Flux_snapshots ...
    ,1-quantiles,3));

figure(99)
z = 0:0.1:20;
J = zeros(5,201);
xJ = zeros(5,201);
testJ = zeros(1,201);
for i = 0:4
    J(i+1,:) = besselj(i,z);
    zJ(i+1,:) = z.*besselj(i,z);
    testJ(1,:) = 0.5*z.*besselj(0,z) + 0.5*z.*besselj(2,z);
end
plot(z,J)
hold on
plot(z,testJ, 'k:', linewidth=3)
%plot(z,zJ)
grid on
legend('J_0','J_1','J_2','J_3','J_4','Location','Best')
title('Bessel Functions of the First Kind for $\nu \in [0, 4]$','interpreter','latex')
xlabel('z','interpreter','latex')
ylabel('$J_\nu(z)$','interpreter','latex')

%% Load coarse-grained data (provided by Cassidy Wagner)

EFlux_CG = ncread("./CoarseGrainedData/multilayerqg_2layer_512_beta_5.0_Ld_0.35_coarsegrain_results.nc", 'EFlux_CG');

PotentialEnstrophyFlux_CG = ncread("./CoarseGrainedData/multilayerqg_2layer_512_beta_5.0_Ld_0.35_coarsegrain_results.nc", 'PotEnstrophy_CG');

EnstrophyFlux_CG = ncread("./CoarseGrainedData/multilayerqg_2layer_512_beta_5.0_Ld_0.35_coarsegrain_results.nc", 'Enstrophy_CG');

K_CG = ncread("./CoarseGrainedData/multilayerqg_2layer_512_beta_5.0_Ld_0.35_coarsegrain_results.nc", 'K_coarse_grain');

K_CG = K_CG/2; % Convert to wavenumber = pi/r rather than 2pi/r (empirically better fit)

alpha = 0.5; % Greyness of coarse-grained flux lines


%% Extract relevant structure functions

Zonal.velveladv = Zonal.uuadv+Zonal.vvadv;
Meridional.velveladv = Meridional.uuadv+Meridional.vvadv;
Diagonal.velveladv = Diagonal.uuadv+Diagonal.vvadv;
OffDiagonal.velveladv = OffDiagonal.uuadv+OffDiagonal.vvadv;

Along_beta = Diagonal;
Across_beta = OffDiagonal;
Intermediate_1 = Zonal;
Intermediate_2 = Meridional;
R_beta = Along_beta.R_DIAGONAL;
R_Intermediate = Along_beta.R;


%% Plot Bessel function estimates of fluxes under isotropic assumption

R_spacing = R_beta(2)-R_beta(1);
R_spacing_Z = R_Intermediate(2)-R_Intermediate(1);

K = Spectral_Flux.Wavenumber;
%K = 2*pi./R_beta; 
%K = J_1_maximum./R_beta;

for layer_no=1:1 % 2 for both layers
    for ii=1:size(K,1)
        % This loop iterates over wavenumber (K) values to find the flux that
        % would be estimated by various structure functions.
        % These estimates utilize Bessel function integrals over all
        % separation distances. Note that first run-through uses diagonal
        % ("Along_beta") separation vectors

        % Enstrophy flux estimates. This example uses the advective
        % vorticity structure function

        % First index represents the upper limit of
        %  separation distance (R) used for the interval. Optimally should
        %  use the maximum R available, but for some purposes it is useful
        %  to know what flux you would estimate from an incomplete
        %  structure function dataset (only measuring up to a specific R)
        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Along_beta.vortvortadv(:,layer_no).*besselj(1,K(ii)*R_beta));

        % This is the enstrophy flux estimate
        bessel_enstrophy_adv_vort_D(ii) = bessel_dummy_integral(end,ii);
        % This is the estimate using a considering all SFs < R(ii)
        % Or equivalently, information from wavenumbers >2*pi/R(ii)
        bessel_enstrophy_adv_vort_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

        % Now repeat for advective q SF estimate of enstrophy flux
        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Along_beta.qqadv(:,layer_no).*besselj(1,K(ii)*R_beta));
        bessel_enstrophy_adv_q_D(ii) = bessel_dummy_integral(end,ii);
        bessel_enstrophy_adv_q_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

        % Now repeat for advective velocity SF estimate of enstrophy flux
        bessel_dummy_integral_D(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing,Along_beta.velveladv(:,layer_no).* (2*besselj(2,K(ii)*R_beta)./(K(ii)*R_beta) - besselj(1,K(ii)*R_beta) ));
        bessel_enstrophy_adv_vel_D(ii) = bessel_dummy_integral(end,ii);
        bessel_enstrophy_adv_vel_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

        % Third-order q-q-Longitudinal velocity
        bessel_dummy_integral_D(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing, ...
            Along_beta.qqL(:,layer_no).*besselj(2,K(ii)*R_beta));
        bessel_enstrophy_Lqq_D(ii) = bessel_dummy_integral(end,ii);
        bessel_enstrophy_Lqq_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

        % Third-order vorticity-vorticity-Longitudinal velocity
        bessel_dummy_integral_D(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing, ...
            Along_beta.vortvortL(:,layer_no).*besselj(2,K(ii)*R_beta));
        bessel_enstrophy_Lvortvort_D(ii) = bessel_dummy_integral(end,ii);
        bessel_enstrophy_Lvortvort_cutoff_D(ii) = bessel_dummy_integral(ii,ii);



        % Now repeat for energy fluxes using:

        % Advective velocity SF
        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Along_beta.velveladv(:,layer_no).*besselj(1,K(ii)*R_beta));
        bessel_energy_adv_vel_D(ii) = bessel_dummy_integral(end,ii);
        bessel_energy_adv_vel_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

        % Third-order longitudinal
        bessel_dummy_integral(:,ii) = -(K(ii)^3)*(1/12)*cumtrapz(R_spacing,Along_beta.LLL(:,layer_no).*besselj(3,K(ii)*R_beta).*R_beta);
        bessel_energy_LLL_D(ii) = bessel_dummy_integral(end,ii);
        bessel_energy_LLL_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

        % Third-order total velocity
        bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing, ...
            (Along_beta.LLL(:,layer_no) + Along_beta.TTL(:,layer_no)).*besselj(2,K(ii)*R_beta));
        bessel_energy_Luu_D(ii) = bessel_dummy_integral(end,ii);
        bessel_energy_Luu_cutoff_D(ii) = bessel_dummy_integral(ii,ii);


        % Repeat all for zonal direction ("Intermediate_1")

        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.vortvortadv(:,layer_no).*besselj(1,K(ii)*R_Intermediate));
        bessel_enstrophy_adv_vort_Z(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.qqadv(:,layer_no).*besselj(1,K(ii)*R_Intermediate));
        bessel_enstrophy_adv_q_Z(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing_Z,Intermediate_1.velveladv(:,layer_no).* ...
            (2*besselj(2,K(ii)*R_Intermediate)./(K(ii)*R_Intermediate) - besselj(1,K(ii)*R_Intermediate) ));
        bessel_enstrophy_adv_vel_Z(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing_Z, ...
            Intermediate_1.qqL(:,layer_no).*besselj(2,K(ii)*R_Intermediate));
        bessel_enstrophy_Lqq_Z(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing_Z, ...
            Intermediate_1.vortvortL(:,layer_no).*besselj(2,K(ii)*R_Intermediate));
        bessel_enstrophy_Lvortvort_Z(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.velveladv(:,layer_no).* ...
            besselj(1,K(ii)*R_Intermediate));
        bessel_energy_adv_vel_Z(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^3)*(1/12)*cumtrapz(R_spacing_Z,Intermediate_1.LLL(:,layer_no).* ...
            besselj(3,K(ii)*R_Intermediate).*R_Intermediate);
        bessel_energy_LLL_Z(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing_Z, ...
            (Intermediate_1.LLL(:,layer_no) + Intermediate_1.TTL(:,layer_no)).*besselj(2,K(ii)*R_Intermediate));
        bessel_energy_Luu_Z(ii) = bessel_dummy_integral(end,ii);

        % Repeat all for meridional direction ("Intermediate_2")

        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.vortvortadv(:,layer_no).*besselj(1,K(ii)*R_Intermediate));
        bessel_enstrophy_adv_vort_M(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.qqadv(:,layer_no).*besselj(1,K(ii)*R_Intermediate));
        bessel_enstrophy_adv_q_M(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing_Z,Intermediate_2.velveladv(:,layer_no).* ...
            (2*besselj(2,K(ii)*R_Intermediate)./(K(ii)*R_Intermediate) - besselj(1,K(ii)*R_Intermediate) ));
        bessel_enstrophy_adv_vel_M(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing_Z, ...
            Intermediate_2.qqL(:,layer_no).*besselj(2,K(ii)*R_Intermediate));
        bessel_enstrophy_Lqq_M(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing_Z, ...
            Intermediate_2.vortvortL(:,layer_no).*besselj(2,K(ii)*R_Intermediate));
        bessel_enstrophy_Lvortvort_M(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.velveladv(:,layer_no).* ...
            besselj(1,K(ii)*R_Intermediate));
        bessel_energy_adv_vel_M(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^3)*(1/12)*cumtrapz(R_spacing_Z,Intermediate_2.LLL(:,layer_no).* ...
            besselj(3,K(ii)*R_Intermediate).*R_Intermediate);
        bessel_energy_LLL_M(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing_Z, ...
            (Intermediate_2.LLL(:,layer_no) + Intermediate_2.TTL(:,layer_no)).*besselj(2,K(ii)*R_Intermediate));
        bessel_energy_Luu_M(ii) = bessel_dummy_integral(end,ii);

        % Repeat all for off-diagonal direction ("Across_beta")

        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.vortvortadv(:,layer_no).*besselj(1,K(ii)*R_beta));
        bessel_enstrophy_adv_vort_OD(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.qqadv(:,layer_no).*besselj(1,K(ii)*R_beta));
        bessel_enstrophy_adv_q_OD(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing,Across_beta.velveladv(:,layer_no).* ...
            (2*besselj(2,K(ii)*R_beta)./(K(ii)*R_beta) - besselj(1,K(ii)*R_beta) ));
        bessel_enstrophy_adv_vel_OD(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing, ...
            Across_beta.qqL(:,layer_no).*besselj(2,K(ii)*R_beta));
        bessel_enstrophy_Lqq_OD(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing, ...
            Across_beta.vortvortL(:,layer_no).*besselj(2,K(ii)*R_beta));
        bessel_enstrophy_Lvortvort_OD(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.velveladv(:,layer_no).* ...
            besselj(1,K(ii)*R_beta));
        bessel_energy_adv_vel_OD(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^3)*(1/12)*cumtrapz(R_spacing,Across_beta.LLL(:,layer_no).* ...
            besselj(3,K(ii)*R_beta).*R_beta);
        bessel_energy_LLL_OD(ii) = bessel_dummy_integral(end,ii);

        bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing, ...
            (Across_beta.LLL(:,layer_no) + Across_beta.TTL(:,layer_no)).*besselj(2,K(ii)*R_beta));
        bessel_energy_Luu_OD(ii) = bessel_dummy_integral(end,ii);

    end



%% QG Enstrophy Flux plots

h7=figure(10+layer_no)
set(h7,'Position',[10 10 2000 1000])

ymin = -1e5;
ymax = 5e5;
if layer_no==2
    ymin = -0.5e4;
    ymax = 1.5e4;
end

% Advective vorticity SFs

subplot(2,3,4)
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux(:,layer_no), 'k-', 'Linewidth', 3);
hold on
Bessel_SFadv_Flux = semilogx(K,bessel_enstrophy_adv_q_Z,'r-', 'Linewidth', 2);
semilogx(K,bessel_enstrophy_adv_q_M,'r-', 'Linewidth', 2);
semilogx(K,bessel_enstrophy_adv_q_D,'r-', 'Linewidth', 2);
semilogx(K,bessel_enstrophy_adv_q_OD,'r-', 'Linewidth', 2);
SFvort_Flux = semilogx(J_1_maximum./R_beta,-0.5*Along_beta.qqadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_beta,-0.5*Across_beta.qqadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_1.qqadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_2.qqadv(:,layer_no),'b-', 'Linewidth', 0.5);
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux(:,layer_no), 'k-', 'Linewidth', 5);
plot([Rossby_deformation_k, Rossby_deformation_k], [ymin, ymax], 'k:', 'Linewidth', 2)
ylabel('\Pi_K^q')
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)
xlabel('Wavenumber K')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
txt = {'SF_{Aq}-based', '\Pi_K^q estimates'};
text(4, 0.8e5 ,txt, fontsize=size_of_font*0.8)
set(gca,'fontsize', size_of_font);

% Advective velocity SFs

subplot(2,3,5)
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.vortvort_Flux(:,layer_no), 'k-', 'Linewidth', 3);
hold on
Bessel_SFadv_Flux = semilogx(K,bessel_enstrophy_adv_vort_Z,'r-', 'Linewidth', 2);
semilogx(K,bessel_enstrophy_adv_vort_M,'r-', 'Linewidth', 2);
semilogx(K,bessel_enstrophy_adv_vort_D,'r-', 'Linewidth', 2);
semilogx(K,bessel_enstrophy_adv_vort_OD,'r-', 'Linewidth', 2);
SFvort_Flux = semilogx(J_1_maximum./R_beta,-0.5*Along_beta.vortvortadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_beta,-0.5*Across_beta.vortvortadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_1.vortvortadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_2.vortvortadv(:,layer_no),'b-', 'Linewidth', 0.5);
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.vortvort_Flux(:,layer_no), 'k-', 'Linewidth', 5);
plot(2*pi./R_beta, ymin, 'k.')
plot([Rossby_deformation_k, Rossby_deformation_k], [ymin, ymax], 'k:', 'Linewidth', 2)
% patch(Patch_k_Bounds, ...
%     [Spectral_Flux.qq_Flux_max(:,layer_no); ...
%     flipud(Spectral_Flux.qq_Flux_min(:,layer_no))]', ...
%     'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)
xlabel('Wavenumber K')
ylabel('\Pi_K^{\omega}')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
txt = {'SF_{A \omega}-based', '\Pi_K^{\omega} estimates'};
text(1, 4e5 ,txt, fontsize=size_of_font*0.8)
set(gca,'fontsize', size_of_font);

subplot(2,3,6)
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.vortvort_Flux(:,layer_no), 'k-', 'Linewidth', 3);
hold on
semilogx(K,bessel_enstrophy_Lvortvort_Z,'r-', 'Linewidth', 2);
semilogx(K,bessel_enstrophy_Lvortvort_M,'r-', 'Linewidth', 2);
semilogx(K,bessel_enstrophy_Lvortvort_D,'r-', 'Linewidth', 2);
semilogx(K,bessel_enstrophy_Lvortvort_OD,'r-', 'Linewidth', 2);
semilogx(J_2_maximum./R_beta,-0.5*Along_beta.vortvortL(:,layer_no)./(R_beta),'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_beta,-0.5*Across_beta.vortvortL(:,layer_no)./(R_beta),'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_Intermediate,-0.5*Intermediate_1.vortvortL(:,layer_no)./(R_Intermediate),'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_Intermediate,-0.5*Intermediate_2.vortvortL(:,layer_no)./(R_Intermediate),'b-', 'Linewidth', 0.5);
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.vortvort_Flux(:,layer_no), 'k-', 'Linewidth', 5);
plot([Rossby_deformation_k, Rossby_deformation_k], [ymin, ymax], 'k:', 'Linewidth', 2)
% patch(Patch_k_Bounds, ...
%     [Spectral_Flux.qq_Flux_max(:,layer_no); ...
%     flipud(Spectral_Flux.qq_Flux_min(:,layer_no))]', ...
%     'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)
xlabel('Wavenumber K')
ylabel('\Pi_K^{\omega}')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
txt = {'SF_{L \omega \omega}-based', '\Pi_K^{\omega} estimates'};
text(1, 4e5 ,txt, fontsize=size_of_font*0.8)
set(gca,'fontsize', size_of_font);

%% QG Energy flux plots
%h7=figure(15+layer_no)
%set(h7,'Position',[10 10 2000 500])

ymin = -600;
ymax = 600;
if layer_no==2
    ymin = -150;
    ymax = 50;
end

% Advective velocity SFs

subplot(2,3,1)
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux(:,layer_no) ,'k-', 'Linewidth', 3);
hold on
Bessel_SFadv_Flux = semilogx(K,bessel_energy_adv_vel_Z,'r-', 'Linewidth', 2);
semilogx(K,bessel_energy_adv_vel_M,'r-', 'Linewidth', 2);
semilogx(K,bessel_energy_adv_vel_D,'r-', 'Linewidth', 2);
semilogx(K,bessel_energy_adv_vel_OD,'r-', 'Linewidth', 2);
SFadv_Flux = semilogx(J_1_maximum./R_beta,-0.5*Along_beta.velveladv(:,layer_no),'b-', 'Linewidth', 0.5);
SFadv_Flux = semilogx(J_1_maximum./R_beta,-0.5*Across_beta.velveladv(:,layer_no),'b-', 'Linewidth', 0.5);
SFadv_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_1.velveladv(:,layer_no),'b-', 'Linewidth', 0.5);
SFadv_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_2.velveladv(:,layer_no),'b-', 'Linewidth', 0.5);
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux(:,layer_no) ,'k-', 'Linewidth', 5);
plot([Rossby_deformation_k, Rossby_deformation_k], [ymin, ymax], 'k:', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)
ylabel('\Pi_K^u')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
legend([Energy_Flux Bessel_SFadv_Flux ...
     SFadv_Flux],'Fourier-estimated Flux','Bessel Method', ...
    'Traditional Method','Location','North');
txt = {'SF_{Au}-based', '\Pi_K^u estimates'};
text(50, -350 ,txt, fontsize=size_of_font*0.8)
set(gca,'fontsize', size_of_font);


% Third-order Lvv and Lvortvort SFs

subplot(2,3,2)
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux(:,layer_no) ,'k-', 'Linewidth', 3);
hold on
semilogx(K,bessel_energy_Luu_Z,'r-', 'Linewidth', 2);
semilogx(K,bessel_energy_Luu_M,'r-', 'Linewidth', 2);
semilogx(K,bessel_energy_Luu_D,'r-', 'Linewidth', 2);
semilogx(K,bessel_energy_Luu_OD,'r-', 'Linewidth', 2);
semilogx(J_2_maximum./R_beta,-0.5*(Along_beta.LLL(:,layer_no) + Along_beta.TTL(:,layer_no))./R_beta,'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_beta,-0.5*(Across_beta.LLL(:,layer_no) + Across_beta.TTL(:,layer_no))./R_beta,'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_Intermediate,-0.5*(Intermediate_1.LLL(:,layer_no) + Intermediate_1.TTL(:,layer_no))./R_Intermediate,'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_Intermediate,-0.5*(Intermediate_2.LLL(:,layer_no) + Intermediate_2.TTL(:,layer_no))./R_Intermediate,'b-', 'Linewidth', 0.5);
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux(:,layer_no) ,'k-', 'Linewidth', 5);
plot([Rossby_deformation_k, Rossby_deformation_k], [ymin, ymax], 'k:', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)
ylim([ymin, ymax]);
xlim([xmin, xmax]);
ylabel('\Pi_K')
txt = {'SF_{Luu}-based', '\Pi_K^u estimates'};
text(50, -350 ,txt, fontsize=size_of_font*0.8)
set(gca,'fontsize', size_of_font);

% Third-order LLL SF

subplot(2,3,3)
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux(:,layer_no) ,'k-', 'Linewidth', 3);
hold on
semilogx(K,bessel_energy_LLL_Z,'r-', 'Linewidth', 2);
semilogx(K,bessel_energy_LLL_M,'r-', 'Linewidth', 2);
semilogx(K,bessel_energy_LLL_D,'r-', 'Linewidth', 2);
semilogx(K,bessel_energy_LLL_OD,'r-', 'Linewidth', 2);
semilogx(J_3_maximum./R_beta,-(2/3)*Along_beta.LLL(:,layer_no)./R_beta,'b-', 'Linewidth', 0.5);
semilogx(J_3_maximum./R_beta,-(2/3)*Across_beta.LLL(:,layer_no)./R_beta,'b-', 'Linewidth', 0.5);
semilogx(J_3_maximum./R_Intermediate,-(2/3)*Intermediate_1.LLL(:,layer_no)./R_Intermediate,'b-', 'Linewidth', 0.5);
semilogx(J_3_maximum./R_Intermediate,-(2/3)*Intermediate_2.LLL(:,layer_no)./R_Intermediate,'b-', 'Linewidth', 0.5);
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux(:,layer_no) ,'k-', 'Linewidth', 5);
plot([Rossby_deformation_k, Rossby_deformation_k], [ymin, ymax], 'k:', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)
ylim([ymin, ymax]);
xlim([xmin, xmax]);
ylabel('\Pi_K')
txt = {'SF_{LLL}-based', '\Pi_K^u estimates'};
text(1, 400 ,txt, fontsize=size_of_font*0.8)
set(gca,'fontsize', size_of_font);

% Create simple plots for main paper text

h7=figure(20+layer_no)
set(h7,'Position',[10 10 2000 500])

ymin = -600;
ymax = 50;
if layer_no==2
    ymin = -150;
    ymax = 50;
end

% Kinetic Energy Plot

bessel_energy_adv_vel_alldirections = ...
    cat(3,bessel_energy_adv_vel_Z, bessel_energy_adv_vel_M, ...
    bessel_energy_adv_vel_D, bessel_energy_adv_vel_OD);

patch_upper_bound = mean(bessel_energy_adv_vel_alldirections, 3) + ...
    std(bessel_energy_adv_vel_alldirections, 1, 3);

patch_lower_bound = mean(bessel_energy_adv_vel_alldirections, 3) - ...
    std(bessel_energy_adv_vel_alldirections, 1, 3);

subplot(1,3,1)
Bessel_SFveladv_Flux = semilogx(K, ...
    mean(bessel_energy_adv_vel_alldirections, 3),'r-', 'Linewidth', 2);
hold on

patch(Patch_k_Bounds, ...
    [patch_upper_bound'; ...
    flipud(patch_lower_bound')]', ...
    'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2)

SFadv_Flux = semilogx(J_1_maximum./R_beta,-0.5*Along_beta.velveladv(:,layer_no),'r--', 'Linewidth', 1);

Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux(:,layer_no) ,'k-', 'Linewidth', 5);

CG_Flux = semilogx(K_CG, EFlux_CG, 'r:', 'Linewidth', 2);

plot([Rossby_deformation_k, Rossby_deformation_k], [ymin, ymax], 'k--', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)
xlabel('Wavenumber K')
ylabel('\Pi_K^u')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
title("Kinetic Energy")
%title("Kinetic Energy (Layer "+layer_no+")")
set(gca,'fontsize', size_of_font);
legend([Energy_Flux Bessel_SFveladv_Flux SFadv_Flux CG_Flux],'Fourier-estimated Flux','Bessel Method: SF_{Au}', ...
    'Traditional Method: SF_{Au}','Coarse-graining ({\bfu\rm})','Location', 'SouthEast', 'fontsize', legend_font_size);


ymin = -2e5;
ymax = 5e5;
if layer_no==2
    ymin = -1e4;
    ymax = 4e4;
end

bessel_enstrophy_adv_q_alldirections = ...
    cat(3,bessel_enstrophy_adv_q_Z, bessel_enstrophy_adv_q_M, ...
    bessel_enstrophy_adv_q_D, bessel_enstrophy_adv_q_OD);

bessel_enstrophy_adv_vort_alldirections = ...
    cat(3,bessel_enstrophy_adv_vort_Z, bessel_enstrophy_adv_vort_M, ...
    bessel_enstrophy_adv_vort_D, bessel_enstrophy_adv_vort_OD);

bessel_enstrophy_adv_vel_alldirections = ...
    cat(3,bessel_enstrophy_adv_vel_Z, bessel_enstrophy_adv_vel_M, ...
    bessel_enstrophy_adv_vel_D, bessel_enstrophy_adv_vel_OD);

% QG Potential Enstrophy Plot

subplot(1,3,2)
Bessel_SFqadv_Flux = semilogx(K, ...
    mean(bessel_enstrophy_adv_q_alldirections, 3),'r-', 'Linewidth', 2);
hold on

patch_upper_bound = mean(bessel_enstrophy_adv_q_alldirections, 3) + ...
    std(bessel_enstrophy_adv_q_alldirections, 1, 3);
patch_lower_bound = mean(bessel_enstrophy_adv_q_alldirections, 3) - ...
    std(bessel_enstrophy_adv_q_alldirections, 1, 3);
patch(Patch_k_Bounds, ...
    [patch_upper_bound'; ...
    flipud(patch_lower_bound')]', ...
    'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2)

SF_qadv_Flux = semilogx(J_1_maximum./R_beta,-0.5*Along_beta.qqadv(:,layer_no),'r--', 'Linewidth', 1);

CG_Flux = semilogx(K_CG, PotentialEnstrophyFlux_CG, 'r:', 'Linewidth', 2);

QG_Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux(:,layer_no), 'k-', 'Linewidth', 5);
plot([Rossby_deformation_k, Rossby_deformation_k], [ymin, ymax], 'k--', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)
ylabel('\Pi^q_K')
legend([QG_Enstrophy_Flux Bessel_SFqadv_Flux SF_qadv_Flux CG_Flux],'Fourier-estimated Flux', 'Bessel Method: SF_{Aq}', ...
    'Traditional Method: SF_{Aq}', 'Coarse-graining (q)', ...
    'Location','South', 'fontsize', legend_font_size);
xlabel('Wavenumber K')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
title("Potential Enstrophy")
%title("Potential Enstrophy (Layer "+layer_no+")")
set(gca,'fontsize', size_of_font);

% Enstrophy Plot

subplot(1,3,3)
Bessel_SFvortadv_Flux = semilogx(K, ...
    mean(bessel_enstrophy_adv_vort_alldirections, 3),'r-', 'Linewidth', 2);
hold on

patch_upper_bound = mean(bessel_enstrophy_adv_vort_alldirections, 3) + ...
    std(bessel_enstrophy_adv_vort_alldirections, 1, 3);
patch_lower_bound = mean(bessel_enstrophy_adv_vort_alldirections, 3) - ...
    std(bessel_enstrophy_adv_vort_alldirections, 1, 3);
patch(Patch_k_Bounds, ...
    [patch_upper_bound'; ...
    flipud(patch_lower_bound')]', ...
    'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2)

Bessel_SFveladv_Flux = semilogx(K, ...
    mean(bessel_enstrophy_adv_vel_alldirections, 3),'b-', 'Linewidth', 2);
patch_upper_bound = mean(bessel_enstrophy_adv_vel_alldirections, 3) + ...
    std(bessel_enstrophy_adv_vel_alldirections, 1, 3);
patch_lower_bound = mean(bessel_enstrophy_adv_vel_alldirections, 3) - ...
    std(bessel_enstrophy_adv_vel_alldirections, 1, 3);
patch(Patch_k_Bounds, ...
    [patch_upper_bound'; ...
    flipud(patch_lower_bound')]', ...
    'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2)

SF_vortadv_Flux = semilogx(J_1_maximum./R_beta,-0.5*Along_beta.vortvortadv(:,layer_no),'r--', 'Linewidth', 1);
SF_veladv_Flux = semilogx(J_1_maximum./R_beta, 2*Along_beta.velveladv(:,layer_no)./(R_beta.^2),'b--', 'Linewidth', 1);

CG_Flux = semilogx(K_CG, EnstrophyFlux_CG, 'r:', 'Linewidth', 2);

Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.vortvort_Flux(:,layer_no), 'k-', 'Linewidth', 5);
plot([Rossby_deformation_k, Rossby_deformation_k], [ymin, ymax], 'k--', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

semilogx(K, mean(bessel_enstrophy_adv_vort_alldirections, 3),'r-', 'Linewidth', 2);
ylabel('\Pi^{\omega}_K')
legend([Enstrophy_Flux Bessel_SFvortadv_Flux Bessel_SFveladv_Flux SF_vortadv_Flux SF_veladv_Flux CG_Flux], 'Fourier-estimated Flux', 'Bessel Method: SF_{A\omega}' , ...
    'Bessel Method: SF_{Au}' , 'Traditional Method: SF_{A\omega}' , 'Traditional Method: SF_{Au}' , 'Coarse-graining (\omega)', ...
    'Location','NorthWest', 'fontsize', legend_font_size);
xlabel('Wavenumber K')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
title("Enstrophy")
%title("Enstrophy (Layer "+layer_no+")")
set(gca,'fontsize', size_of_font);



end