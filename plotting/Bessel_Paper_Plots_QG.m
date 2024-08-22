%% Bessel function tests

clear all

size_of_font = 20;
J_1_maximum = 2;
J_2_maximum = 3.1;
J_3_maximum = 4.3;

experiment_name = "multilayerqg_2layer_1024"; %["SB","WB","SBd"];

quantiles = 0.75;

no_of_exps = size(experiment_name,2);

%k_Energy_Cascade_End = 0.9e1;
%r_Energy_Cascade_End = 2*pi/k_Energy_Cascade_End;

k_Enstrophy_Cascade_End = 4e2;
r_Enstrophy_Cascade_End = 2*pi/k_Enstrophy_Cascade_End;

k_forcing_estimate = 1e1;
r_forcing_estimate = 2*pi/k_forcing_estimate;
path = "../../Paper_Anisotropic_2D_Structure_Functions/analysis/";

Diagonal = load('Structure_Functions_'+experiment_name+'_Diagonal.mat');
Meridional = load('Structure_Functions_'+experiment_name+'_Meridional.mat');
Zonal = load('Structure_Functions_'+experiment_name+'_Zonal.mat');
OffDiagonal = load('Structure_Functions_'+experiment_name+'_Off-Diagonal.mat');


%% Plot spectral fluxes

Spectral_Flux = load('Spectral_Fluxes_'+experiment_name+'.mat');

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
for layer_no=1:2
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

h7=figure(11+layer_no)
set(h7,'Position',[10 10 2000 500])

ymin = -10;
    ymax = 20;
    if layer_no==2 
        ymin = -0.05;
        ymax = 0.15;
    end

% Advective vorticity SFs

subplot(1,3,1)
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux(:,layer_no), 'k-', 'Linewidth', 3);
hold on
Bessel_SFadv_Flux = semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_q_Z,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_q_M,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_q_D,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_q_OD,'r-', 'Linewidth', 2);
SFvort_Flux = semilogx(J_1_maximum./R_beta,-0.5*Along_beta.qqadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_beta,-0.5*Across_beta.qqadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_1.qqadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_2.qqadv(:,layer_no),'b-', 'Linewidth', 0.5);
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux(:,layer_no), 'k-', 'Linewidth', 5);
patch(Patch_k_Bounds, ...
    [Spectral_Flux.qq_Flux_max(:,layer_no); ...
    flipud(Spectral_Flux.qq_Flux_min(:,layer_no))]', ...
    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
ylabel('Spectral qq/2 flux (s^{-3})')
legend([Enstrophy_Flux Bessel_SFadv_Flux ...
    SFvort_Flux],'qq/2 Flux','Bessel Method', ...
    'Basic Method' ,'Location','NorthWest');
xlabel('Wavenumber K (rad m^{-1})')
ylim([ymin, ymax]);
title('SF_{Aq}')
set(gca,'fontsize', size_of_font);

% Advective velocity SFs

subplot(1,3,2)
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux(:,layer_no), 'k-', 'Linewidth', 3);
hold on
Bessel_SFadv_Flux = semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vort_Z,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vort_M,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vort_D,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_adv_vort_OD,'r-', 'Linewidth', 2);
SFvort_Flux = semilogx(J_1_maximum./R_beta,-0.5*Along_beta.vortvortadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_beta,-0.5*Across_beta.vortvortadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_1.vortvortadv(:,layer_no),'b-', 'Linewidth', 0.5);
SFvort_Flux = semilogx(J_1_maximum./R_Intermediate,-0.5*Intermediate_2.vortvortadv(:,layer_no),'b-', 'Linewidth', 0.5);
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux(:,layer_no), 'k-', 'Linewidth', 5);
plot(2*pi./R_beta, ymin, 'k.')
patch(Patch_k_Bounds, ...
    [Spectral_Flux.qq_Flux_max(:,1); ...
    flipud(Spectral_Flux.qq_Flux_min(:,1))]', ...
    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
xlabel('Wavenumber K (rad m^{-1})')
ylim([ymin, ymax]);
title('SF_{A\omega}')
set(gca,'fontsize', size_of_font);

subplot(1,3,3)
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux(:,layer_no), 'k-', 'Linewidth', 3);
hold on
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_Lvortvort_Z,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_Lvortvort_M,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_Lvortvort_D,'r-', 'Linewidth', 2);
semilogx(Spectral_Flux.Wavenumber,bessel_enstrophy_Lvortvort_OD,'r-', 'Linewidth', 2);
semilogx(J_2_maximum./R_beta,-0.5*Along_beta.vortvortL(:,layer_no)./(R_beta),'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_beta,-0.5*Across_beta.vortvortL(:,layer_no)./(R_beta),'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_Intermediate,-0.5*Intermediate_1.vortvortL(:,layer_no)./(R_Intermediate),'b-', 'Linewidth', 0.5);
semilogx(J_2_maximum./R_Intermediate,-0.5*Intermediate_2.vortvortL(:,layer_no)./(R_Intermediate),'b-', 'Linewidth', 0.5);
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux(:,layer_no), 'k-', 'Linewidth', 5);
patch(Patch_k_Bounds, ...
    [Spectral_Flux.qq_Flux_max(:,1); ...
    flipud(Spectral_Flux.qq_Flux_min(:,1))]', ...
    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
xlabel('Wavenumber K (rad m^{-1})')
ylim([ymin, ymax]);
title('SF_{L \omega \omega}')
set(gca,'fontsize', size_of_font);

end