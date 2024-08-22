%% Bessel function paper plots (2D turbulence)

clear all

size_of_font = 20;
J_1_maximum = 2;
J_2_maximum = 3.1;
J_3_maximum = 4.3;
J_2_3_function_minimum = 1.5;

experiment_name = "StrongBeta"; %["SB","WB","SBd"];

quantiles = 0.75;

no_of_exps = size(experiment_name,2);

%k_Energy_Cascade_End = 0.9e1;
%r_Energy_Cascade_End = 2*pi/k_Energy_Cascade_End;

k_Enstrophy_Cascade_End = 4e2;
r_Enstrophy_Cascade_End = 2*pi/k_Enstrophy_Cascade_End;

k_forcing_estimate = 1e1;
r_forcing_estimate = 2*pi/k_forcing_estimate;
path = "../analysis/processed_data/";

Diagonal = load(path+'Structure_Functions_2D_'+experiment_name+'_Diagonal.mat');
Meridional = load(path+'Structure_Functions_2D_'+experiment_name+'_Meridional.mat');
Zonal = load(path+'Structure_Functions_2D_'+experiment_name+'_Zonal.mat');
OffDiagonal = load(path+'Structure_Functions_2D_'+experiment_name+'_Off-Diagonal.mat');

%% Plot spectral fluxes

Spectral_Flux = load(path+'Spectral_Fluxes_2D_'+experiment_name+'.mat');

Patch_k_Bounds = [Spectral_Flux.Wavenumber' fliplr(Spectral_Flux.Wavenumber')];
Patch_k_Bounds_Enstrophy = [k_forcing_estimate, max(Spectral_Flux.Wavenumber), ...
    max(Spectral_Flux.Wavenumber), k_forcing_estimate];
%Patch_k_Bounds_Energy = [k_Energy_Cascade_End, k_forcing_estimate, ...
%    k_forcing_estimate, k_Energy_Cascade_End];

Spectral_Flux.Enstrophy_Flux_max = (quantile(Spectral_Flux.Enstrophy_Flux_snapshots ...
    ,quantiles,3))-Spectral_Flux.Enstrophy_Flux(1);
Spectral_Flux.Enstrophy_Flux_min = (quantile(Spectral_Flux.Enstrophy_Flux_snapshots ...
    ,1-quantiles,3))-Spectral_Flux.Enstrophy_Flux(1);
Spectral_Flux.Energy_Flux_max = (quantile(Spectral_Flux.Energy_Flux_snapshots ...
    ,quantiles,3));
Spectral_Flux.Energy_Flux_min = (quantile(Spectral_Flux.Energy_Flux_snapshots ...
    ,1-quantiles,3));

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
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Along_beta.vortvortadv(:,1).*besselj(1,K(ii)*R_beta));

    % This is the enstrophy flux estimate
    bessel_enstrophy_adv_vort_D(ii) = bessel_dummy_integral(end,ii);
    % This is the estimate using a considering all SFs < R(ii)
    % Or equivalently, information from wavenumbers >2*pi/R(ii)
    bessel_enstrophy_adv_vort_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

    % Now repeat for advective velocity SF estimate of enstrophy flux
    bessel_dummy_integral_D(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing,Along_beta.velveladv(:,1).* (2*besselj(2,K(ii)*R_beta)./(K(ii)*R_beta) - besselj(1,K(ii)*R_beta) ));
    bessel_enstrophy_adv_vel_D(ii) = bessel_dummy_integral(end,ii);
    bessel_enstrophy_adv_vel_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

    % Third-order velocity-velocity-Longitudinal velocity
    bessel_dummy_integral_D(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing, ...
        Along_beta.vortvortL.*besselj(2,K(ii)*R_beta));
    bessel_enstrophy_Lvortvort_D(ii) = bessel_dummy_integral(end,ii);
    bessel_enstrophy_Lvortvort_cutoff_D(ii) = bessel_dummy_integral(ii,ii);



    % Now repeat for energy fluxes using:

    % Advective velocity SF
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Along_beta.velveladv(:,1).*besselj(1,K(ii)*R_beta));
    bessel_energy_adv_vel_D(ii) = bessel_dummy_integral(end,ii);
    bessel_energy_adv_vel_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

    % Third-order longitudinal
    bessel_dummy_integral(:,ii) = -(K(ii)^3)*(1/12)*cumtrapz(R_spacing,Along_beta.LLL(:,1).*besselj(3,K(ii)*R_beta).*R_beta);
    bessel_energy_LLL_D(ii) = bessel_dummy_integral(end,ii);
    bessel_energy_LLL_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

    % Third-order total velocity
    bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing, ...
        (Along_beta.LLL(:,1) + Along_beta.TTL(:,1)).*besselj(2,K(ii)*R_beta));
    bessel_energy_Luu_D(ii) = bessel_dummy_integral(end,ii);
    bessel_energy_Luu_cutoff_D(ii) = bessel_dummy_integral(ii,ii);

    
    % Repeat all for zonal direction ("Intermediate_1")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.vortvortadv(:,1).*besselj(1,K(ii)*R_Intermediate));
    bessel_enstrophy_adv_vort_Z(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing_Z,Intermediate_1.velveladv(:,1).* ...
        (2*besselj(2,K(ii)*R_Intermediate)./(K(ii)*R_Intermediate) - besselj(1,K(ii)*R_Intermediate) ));
    bessel_enstrophy_adv_vel_Z(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing_Z, ...
        Intermediate_1.vortvortL.*besselj(2,K(ii)*R_Intermediate));
    bessel_enstrophy_Lvortvort_Z(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.velveladv(:,1).* ...
        besselj(1,K(ii)*R_Intermediate));
    bessel_energy_adv_vel_Z(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = -(K(ii)^3)*(1/12)*cumtrapz(R_spacing_Z,Intermediate_1.LLL(:,1).* ...
        besselj(3,K(ii)*R_Intermediate).*R_Intermediate);
    bessel_energy_LLL_Z(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing_Z, ...
        (Intermediate_1.LLL(:,1) + Intermediate_1.TTL(:,1)).*besselj(2,K(ii)*R_Intermediate));
    bessel_energy_Luu_Z(ii) = bessel_dummy_integral(end,ii);

    % Repeat all for meridional direction ("Intermediate_2")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.vortvortadv(:,1).*besselj(1,K(ii)*R_Intermediate));
    bessel_enstrophy_adv_vort_M(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing_Z,Intermediate_2.velveladv(:,1).* ...
        (2*besselj(2,K(ii)*R_Intermediate)./(K(ii)*R_Intermediate) - besselj(1,K(ii)*R_Intermediate) ));
    bessel_enstrophy_adv_vel_M(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing_Z, ...
        Intermediate_2.vortvortL.*besselj(2,K(ii)*R_Intermediate));
    bessel_enstrophy_Lvortvort_M(ii) = bessel_dummy_integral(end,ii);
    
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.velveladv(:,1).* ...
        besselj(1,K(ii)*R_Intermediate));
    bessel_energy_adv_vel_M(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = -(K(ii)^3)*(1/12)*cumtrapz(R_spacing_Z,Intermediate_2.LLL(:,1).* ...
        besselj(3,K(ii)*R_Intermediate).*R_Intermediate);
    bessel_energy_LLL_M(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing_Z, ...
        (Intermediate_2.LLL(:,1) + Intermediate_2.TTL(:,1)).*besselj(2,K(ii)*R_Intermediate));
    bessel_energy_Luu_M(ii) = bessel_dummy_integral(end,ii);

    % Repeat all for off-diagonal direction ("Across_beta")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.vortvortadv(:,1).*besselj(1,K(ii)*R_beta));
    bessel_enstrophy_adv_vort_OD(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing,Across_beta.velveladv(:,1).* ...
        (2*besselj(2,K(ii)*R_beta)./(K(ii)*R_beta) - besselj(1,K(ii)*R_beta) ));
    bessel_enstrophy_adv_vel_OD(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing, ...
        Across_beta.vortvortL.*besselj(2,K(ii)*R_beta));
    bessel_enstrophy_Lvortvort_OD(ii) = bessel_dummy_integral(end,ii);
    
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.velveladv(:,1).* ...
        besselj(1,K(ii)*R_beta));
    bessel_energy_adv_vel_OD(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = -(K(ii)^3)*(1/12)*cumtrapz(R_spacing,Across_beta.LLL(:,1).* ...
        besselj(3,K(ii)*R_beta).*R_beta);
    bessel_energy_LLL_OD(ii) = bessel_dummy_integral(end,ii);

    bessel_dummy_integral(:,ii) = -(K(ii)^2)*(1/4)*cumtrapz(R_spacing, ...
        (Across_beta.LLL(:,1) + Across_beta.TTL(:,1)).*besselj(2,K(ii)*R_beta));
    bessel_energy_Luu_OD(ii) = bessel_dummy_integral(end,ii);
    
end

%% 2D Energy flux plots
h7=figure(10)
set(h7,'Position',[10 10 2000 500])

% Advective velocity SFs

subplot(1,3,1)
ymin = -1.5e-5;
ymax = 1e-5;
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Energy_Flux(:,1) ,'k-', 'Linewidth', 3);
hold on
Bessel_SFadv_Flux = semilogx(Spectral_Flux.Wavenumber,bessel_energy_adv_vel_Z,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_energy_adv_vel_M,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_energy_adv_vel_D,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_energy_adv_vel_OD,'r-', 'Linewidth', 2);
SFadv_Flux = semilogx(J_1_maximum./R_beta,-0.5*Along_beta.velveladv(:,1),'b-', 'Linewidth', 0.5);
SFadv_Flux = semilogx(J_1_maximum./R_beta,-0.5*Across_beta.velveladv(:,1),'b-', 'Linewidth', 0.5);
SFadv_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_1.velveladv(:,1),'b-', 'Linewidth', 0.5);
SFadv_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_2.velveladv(:,1),'b-', 'Linewidth', 0.5);
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Energy_Flux(:,1) ,'k-', 'Linewidth', 5);
patch(Patch_k_Bounds, ...
    [Spectral_Flux.Energy_Flux_max(:,1); ...
    flipud(Spectral_Flux.Energy_Flux_min(:,1))]', ...
    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
xlabel('Wavenumber K (rad m^{-1})')
ylabel('Spectral Energy flux (m^2 s^{-3})')
ylim([ymin, ymax]);
legend([Energy_Flux Bessel_SFadv_Flux ...
     SFadv_Flux],'Energy Flux','Bessel Method', ...
    'Basic method','Location','NorthEast');
title('SF_{Au}')
set(gca,'fontsize', size_of_font);


% Third-order Lvv and Lvortvort SFs

subplot(1,3,2)
ymin = -1.5e-5;
ymax = 1e-5;
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Energy_Flux(:,1) ,'k-', 'Linewidth', 3);
hold on
semilogx(Spectral_Flux.Wavenumber,bessel_energy_Luu_Z,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_energy_Luu_M,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_energy_Luu_D,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_energy_Luu_OD,'r-', 'Linewidth', 2);
semilogx(J_2_maximum./R_beta,-0.5*(Along_beta.LLL(:,1) + Along_beta.TTL(:,1))./R_beta,'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_beta,-0.5*(Across_beta.LLL(:,1) + Across_beta.TTL(:,1))./R_beta,'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_Intermediate,-0.5*(Intermediate_1.LLL(:,1) + Intermediate_1.TTL(:,1))./R_Intermediate,'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_Intermediate,-0.5*(Intermediate_2.LLL(:,1) + Intermediate_2.TTL(:,1))./R_Intermediate,'b-', 'Linewidth', 0.5);
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Energy_Flux(:,1) ,'k-', 'Linewidth', 5);
patch(Patch_k_Bounds, ...
    [Spectral_Flux.Energy_Flux_max(:,1); ...
    flipud(Spectral_Flux.Energy_Flux_min(:,1))]', ...
    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
xlabel('Wavenumber K (rad m^{-1})')
ylim([ymin, ymax]);
title('SF_{Luu}')
set(gca,'fontsize', size_of_font);

% Third-order LLL SF

subplot(1,3,3)
ymin = -1.5e-5;
ymax = 1e-5;
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Energy_Flux(:,1) ,'k-', 'Linewidth', 3);
hold on
semilogx(Spectral_Flux.Wavenumber,bessel_energy_LLL_Z,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_energy_LLL_M,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_energy_LLL_D,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_energy_LLL_OD,'r-', 'Linewidth', 2);
semilogx(J_3_maximum./R_beta,-(2/3)*Along_beta.LLL(:,1)./R_beta,'b-', 'Linewidth', 0.5);
semilogx(J_3_maximum./R_beta,-(2/3)*Across_beta.LLL(:,1)./R_beta,'b-', 'Linewidth', 0.5);
semilogx(J_3_maximum./R_Intermediate,-(2/3)*Intermediate_1.LLL(:,1)./R_Intermediate,'b-', 'Linewidth', 0.5);
semilogx(J_3_maximum./R_Intermediate,-(2/3)*Intermediate_2.LLL(:,1)./R_Intermediate,'b-', 'Linewidth', 0.5);
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Energy_Flux(:,1) ,'k-', 'Linewidth', 5);
patch(Patch_k_Bounds, ...
    [Spectral_Flux.Energy_Flux_max(:,1); ...
    flipud(Spectral_Flux.Energy_Flux_min(:,1))]', ...
    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
xlabel('Wavenumber K (rad m^{-1})')
ylim([ymin, ymax]);
title('SF_{LLL}')
set(gca,'fontsize', size_of_font);


%% 2D Enstrophy Flux plots

h7=figure(11)
set(h7,'Position',[10 10 2000 500])

% Advective vorticity SFs

subplot(1,3,1)
ymin = -4e-2;
ymax = 0.12;
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Enstrophy_Flux(:,1), 'k-', 'Linewidth', 3);
hold on
Bessel_SFadv_Flux = semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vort_Z,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vort_M,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vort_D,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vort_OD,'r-', 'Linewidth', 2);
SFvort_Flux = semilogx(J_1_maximum./R_beta,-0.5*Along_beta.vortvortadv(:,1),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_beta,-0.5*Across_beta.vortvortadv(:,1),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_1.vortvortadv(:,1),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_2.vortvortadv(:,1),'b-', 'Linewidth', 0.5);
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Enstrophy_Flux(:,1), 'k-', 'Linewidth', 5);
patch(Patch_k_Bounds, ...
    [Spectral_Flux.Enstrophy_Flux_max(:,1); ...
    flipud(Spectral_Flux.Enstrophy_Flux_min(:,1))]', ...
    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
ylabel('Spectral Enstrophy flux (s^{-3})')
legend([Enstrophy_Flux Bessel_SFadv_Flux ...
    SFvort_Flux],'Enstrophy Flux','Bessel Method', ...
    'Basic Method' ,'Location','NorthWest');
xlabel('Wavenumber K (rad m^{-1})')
ylim([ymin, ymax]);
title('SF_{A\omega}')
set(gca,'fontsize', size_of_font);

% Advective velocity SFs

subplot(1,3,2)
ymin = -4e-2;
ymax = 0.12;
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Enstrophy_Flux(:,1), 'k-', 'Linewidth', 3);
hold on
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vel_Z,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vel_M,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vel_D,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vel_OD,'r-', 'Linewidth', 2);
semilogx(J_2_3_function_minimum./R_beta,2*Along_beta.velveladv(:,1)./(R_beta.^2),'b-', 'Linewidth', 0.5);
semilogx(J_2_3_function_minimum./R_beta,2*Across_beta.velveladv(:,1)./(R_beta.^2),'b-', 'Linewidth', 0.5);
semilogx(J_2_3_function_minimum./R_Intermediate,2*Intermediate_1.velveladv(:,1)./(R_Intermediate.^2),'b-', 'Linewidth', 0.5);
semilogx(J_2_3_function_minimum./R_Intermediate,2*Intermediate_2.velveladv(:,1)./(R_Intermediate.^2),'b-', 'Linewidth', 0.5);
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Enstrophy_Flux(:,1), 'k-', 'Linewidth', 5);
plot(2*pi./R_beta, ymin, 'k.')
patch(Patch_k_Bounds, ...
    [Spectral_Flux.Enstrophy_Flux_max(:,1); ...
    flipud(Spectral_Flux.Enstrophy_Flux_min(:,1))]', ...
    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
xlabel('Wavenumber K (rad m^{-1})')
ylim([ymin, ymax]);
title('SF_{Au}')
set(gca,'fontsize', size_of_font);

subplot(1,3,3)
ymin = -4e-2;
ymax = 0.12;
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Enstrophy_Flux(:,1), 'k-', 'Linewidth', 3);
hold on
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_Lvortvort_Z,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_Lvortvort_M,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_Lvortvort_D,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_Lvortvort_OD,'r-', 'Linewidth', 2);
semilogx(J_2_maximum./R_beta,-0.5*Along_beta.vortvortL(:,1)./(R_beta),'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_beta,-0.5*Across_beta.vortvortL(:,1)./(R_beta),'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_Intermediate,-0.5*Intermediate_1.vortvortL(:,1)./(R_Intermediate),'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_Intermediate,-0.5*Intermediate_2.vortvortL(:,1)./(R_Intermediate),'b-', 'Linewidth', 0.5);
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Enstrophy_Flux(:,1), 'k-', 'Linewidth', 5);
patch(Patch_k_Bounds, ...
    [Spectral_Flux.Enstrophy_Flux_max(:,1); ...
    flipud(Spectral_Flux.Enstrophy_Flux_min(:,1))]', ...
    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
xlabel('Wavenumber K (rad m^{-1})')
ylim([ymin, ymax]);
title('SF_{L \omega \omega}')
set(gca,'fontsize', size_of_font);
