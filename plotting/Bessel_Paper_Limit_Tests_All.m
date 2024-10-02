%% Bessel function paper plots (Tests showing impacts of limited obs. range)

clear all

size_of_font = 14;
J_1_maximum = 2;
J_2_maximum = 3.1;
J_3_maximum = 4.3;
J_2_3_function_minimum = 1.5;

%% Start with 2D turbulence simulation

experiment_name = "StrongBeta"; %["SB","WB","SBd"];

quantiles = 0.75;

no_of_exps = size(experiment_name,2);

%k_Energy_Cascade_End = 0.9e1;
%r_Energy_Cascade_End = 2*pi/k_Energy_Cascade_End;

k_Enstrophy_Cascade_End = 4e2;
r_Enstrophy_Cascade_End = 2*pi/k_Enstrophy_Cascade_End;

k_forcing_estimate = 1e1;
forcing_k = 1e2;
max_energy_cascade_k = 1e1;
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

n_min=[2; 4; 8; 16];
n_max = (1./n_min)*size(R_beta,1)/2;

n_min_D = ceil(n_min/sqrt(2));
n_max_D = floor(n_max/sqrt(2));


%% Calculate 2D turbulence Bessel function estimates of fluxes

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

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vort_rmin_cutoff_D(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_enstrophy_adv_vort_rmax_cutoff_D(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end

    % Now repeat for advective velocity SF estimate of enstrophy flux
    bessel_dummy_integral_D(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing,Along_beta.velveladv(:,1).* (2*besselj(2,K(ii)*R_beta)./(K(ii)*R_beta) - besselj(1,K(ii)*R_beta) ));
    bessel_enstrophy_adv_vel_D(ii) = bessel_dummy_integral(end,ii);
    
    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vel_rmin_cutoff_D(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_enstrophy_adv_vel_rmax_cutoff_D(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end
    


    % Now repeat for energy fluxes using:

    % Advective velocity SF
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Along_beta.velveladv(:,1).*besselj(1,K(ii)*R_beta));
    bessel_energy_adv_vel_D(ii) = bessel_dummy_integral(end,ii);
    %bessel_energy_adv_vel_cutoff_D(ii) = bessel_dummy_integral(ii,ii);
    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_D(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_D(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
        r_min_D(k) = R_beta(n_min_D(k));
        r_max_D(k) = R_beta(n_max_D(k));
    end

    
    % Repeat all for zonal direction ("Intermediate_1")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.vortvortadv(:,1).*besselj(1,K(ii)*R_Intermediate));
    bessel_enstrophy_adv_vort_Z(ii) = bessel_dummy_integral(end,ii);

     for k=1:size(n_min,1)
        bessel_enstrophy_adv_vort_rmin_cutoff_Z(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_enstrophy_adv_vort_rmax_cutoff_Z(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end

    bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing_Z,Intermediate_1.velveladv(:,1).* ...
        (2*besselj(2,K(ii)*R_Intermediate)./(K(ii)*R_Intermediate) - besselj(1,K(ii)*R_Intermediate) ));
    bessel_enstrophy_adv_vel_Z(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vel_rmin_cutoff_Z(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_enstrophy_adv_vel_rmax_cutoff_Z(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.velveladv(:,1).* ...
        besselj(1,K(ii)*R_Intermediate));
    bessel_energy_adv_vel_Z(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_Z(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_Z(k,ii) = bessel_dummy_integral(n_max(k),ii);
        r_min(k) = R_Intermediate(n_min(k));
        r_max(k) = R_Intermediate(n_max(k));
    end


    % Repeat all for meridional direction ("Intermediate_2")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.vortvortadv(:,1).*besselj(1,K(ii)*R_Intermediate));
    bessel_enstrophy_adv_vort_M(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vort_rmin_cutoff_M(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_enstrophy_adv_vort_rmax_cutoff_M(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end

    bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing_Z,Intermediate_2.velveladv(:,1).* ...
        (2*besselj(2,K(ii)*R_Intermediate)./(K(ii)*R_Intermediate) - besselj(1,K(ii)*R_Intermediate) ));
    bessel_enstrophy_adv_vel_M(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vel_rmin_cutoff_M(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_enstrophy_adv_vel_rmax_cutoff_M(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end
    
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.velveladv(:,1).* ...
        besselj(1,K(ii)*R_Intermediate));
    bessel_energy_adv_vel_M(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_M(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_M(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end



    % Repeat all for off-diagonal direction ("Across_beta")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.vortvortadv(:,1).*besselj(1,K(ii)*R_beta));
    bessel_enstrophy_adv_vort_OD(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vort_rmin_cutoff_OD(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_enstrophy_adv_vort_rmax_cutoff_OD(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end

    bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing,Across_beta.velveladv(:,1).* ...
        (2*besselj(2,K(ii)*R_beta)./(K(ii)*R_beta) - besselj(1,K(ii)*R_beta) ));
    bessel_enstrophy_adv_vel_OD(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vel_rmin_cutoff_OD(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_enstrophy_adv_vel_rmax_cutoff_OD(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end
    
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.velveladv(:,1).* ...
        besselj(1,K(ii)*R_beta));
    bessel_energy_adv_vel_OD(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_OD(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_OD(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end
    
end




%% Create simple plots for main paper text

%% Create Energy cascade plots

h7=figure(20)
set(h7,'Position',[10 10 1600 1400])
ymin = -1.5e-5;
ymax = 1e-5;
xmin = 0.7;
xmax = 1000;

% Advective velocity SFs

K_max = J_1_maximum./r_min;
K_min = J_1_maximum./r_max;
for k=1:size(n_min,1)
    K_max_index(k) = size(K(K<K_max(k)),1);
    K_min_index(k) = size(K(K<K_min(k)),1); % Add 1 since first index used is the one after this limit
end

bessel_energy_adv_vel_alldirections = ...
    cat(3,bessel_energy_adv_vel_Z, bessel_energy_adv_vel_M, ...
    bessel_energy_adv_vel_D, bessel_energy_adv_vel_OD);

bessel_energy_adv_vel_alldirections_rmin_cutoff = ...
    cat(3,bessel_energy_adv_vel_rmin_cutoff_Z, bessel_energy_adv_vel_rmin_cutoff_M, ...
    bessel_energy_adv_vel_rmin_cutoff_D, bessel_energy_adv_vel_rmin_cutoff_OD);

bessel_energy_adv_vel_alldirections_rmax_cutoff = ...
    cat(3,bessel_energy_adv_vel_rmax_cutoff_Z, bessel_energy_adv_vel_rmax_cutoff_M, ...
    bessel_energy_adv_vel_rmax_cutoff_D, bessel_energy_adv_vel_rmax_cutoff_OD);


subplot(2,3,1)
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Energy_Flux ,'k-', 'Linewidth', 5);
hold on

smallest_rmin = semilogx(K(1:K_max_index(1)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(1,1:K_max_index(1),:), 3), 'c-', 'Linewidth', 2);
semilogx(K(1:K_max_index(2)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(2,1:K_max_index(2),:), 3), 'b-', 'Linewidth', 2);
semilogx(K(1:K_max_index(3)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(3,1:K_max_index(3),:), 3), 'm-', 'Linewidth', 2);
largest_rmin = semilogx(K(1:K_max_index(4)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(4,1:K_max_index(4),:), 3), 'r-', 'Linewidth', 2);

plot([K_max(1), K_max(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_max(2), K_max(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_max(3), K_max(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_max(4), K_max(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

%semilogx(J_1_maximum./R_beta,-0.5*Along_beta.velveladv(:,1),'r--', 'Linewidth', 1);

%plot([forcing_k, forcing_k], [ymin, ymax], 'k--', 'Linewidth', 2)
%plot([max_energy_cascade_k, max_energy_cascade_k], [ymin, ymax], 'k--', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

%semilogx(K, ...
 %   mean(bessel_energy_adv_vel_alldirections, 3),'r-', 'Linewidth', 2);
ylabel('\Pi_K^u')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
legend([Energy_Flux smallest_rmin largest_rmin], ...
    'KE Flux','SF_{Au}: Smallest r_{min}', 'SF_{Au}: Largest r_{min}', ...
    'Location','NorthWest');
title("2D Turbulence")
set(gca,'fontsize', size_of_font);


h8=figure(21)
set(h8,'Position',[10 10 1400 1400])

subplot(2,3,1)
Energy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Energy_Flux ,'k-', 'Linewidth', 5);
%Bessel_SFveladv_Flux = semilogx(K, ...
%    mean(bessel_energy_adv_vel_D, 3),'r-', 'Linewidth', 2);
hold on

largest_rmax = semilogx(K(K_min_index(1):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(1,K_min_index(1):end,:), 3), 'c-', 'Linewidth', 2);
semilogx(K(K_min_index(2):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(2,K_min_index(2):end,:), 3), 'b-', 'Linewidth', 2);
semilogx(K(K_min_index(3):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(3,K_min_index(3):end,:), 3), 'm-', 'Linewidth', 2);
smallest_rmax = semilogx(K(K_min_index(4):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(4,K_min_index(4):end,:), 3), 'r-', 'Linewidth', 2);

plot([K_min(1), K_min(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_min(2), K_min(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_min(3), K_min(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_min(4), K_min(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

%semilogx(J_1_maximum./R_beta,-0.5*Along_beta.velveladv(:,1),'r--', 'Linewidth', 1);

%plot([forcing_k, forcing_k], [ymin, ymax], 'k--', 'Linewidth', 2)
%plot([max_energy_cascade_k, max_energy_cascade_k], [ymin, ymax], 'k--', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

%semilogx(K, ...
 %   mean(bessel_energy_adv_vel_alldirections, 3),'r-', 'Linewidth', 2);
ylabel('\Pi_K^u')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
legend([Energy_Flux largest_rmax smallest_rmax], ...
    'KE Flux','SF_{Au}: Largest r_{max}', 'SF_{Au}: Smallest r_{max}', ...
    'Location','NorthEast');
title("2D Turbulence")
set(gca,'fontsize', size_of_font);

%% Create enstrophy cascade plots

figure(20)

ymin = -4e-2;
ymax = 0.12;

% Advective vorticity SFs

bessel_enstrophy_adv_vel_alldirections = ...
    cat(3,bessel_enstrophy_adv_vel_Z, bessel_enstrophy_adv_vel_M, ...
    bessel_enstrophy_adv_vel_D, bessel_enstrophy_adv_vel_OD);

bessel_enstrophy_adv_vel_alldirections_rmin_cutoff = ...
    cat(3,bessel_enstrophy_adv_vel_rmin_cutoff_Z, bessel_enstrophy_adv_vel_rmin_cutoff_M, ...
    bessel_enstrophy_adv_vel_rmin_cutoff_D, bessel_enstrophy_adv_vel_rmin_cutoff_OD);

bessel_enstrophy_adv_vel_alldirections_rmax_cutoff = ...
    cat(3,bessel_enstrophy_adv_vel_rmax_cutoff_Z, bessel_enstrophy_adv_vel_rmax_cutoff_M, ...
    bessel_enstrophy_adv_vel_rmax_cutoff_D, bessel_enstrophy_adv_vel_rmax_cutoff_OD);

bessel_enstrophy_adv_vort_alldirections_rmin_cutoff = ...
    cat(3,bessel_enstrophy_adv_vort_rmin_cutoff_Z, bessel_enstrophy_adv_vort_rmin_cutoff_M, ...
    bessel_enstrophy_adv_vort_rmin_cutoff_D, bessel_enstrophy_adv_vort_rmin_cutoff_OD);

bessel_enstrophy_adv_vort_alldirections_rmax_cutoff = ...
    cat(3,bessel_enstrophy_adv_vort_rmax_cutoff_Z, bessel_enstrophy_adv_vort_rmax_cutoff_M, ...
    bessel_enstrophy_adv_vort_rmax_cutoff_D, bessel_enstrophy_adv_vort_rmax_cutoff_OD);

subplot(2,3,4)
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Enstrophy_Flux , 'k-', 'Linewidth', 5);
hold on

au = semilogx(K(1:K_max_index(1)), mean(bessel_enstrophy_adv_vel_alldirections_rmin_cutoff(1,1:K_max_index(1),:), 3), 'c-', 'Linewidth', 2);
semilogx(K(1:K_max_index(2)), mean(bessel_enstrophy_adv_vel_alldirections_rmin_cutoff(2,1:K_max_index(2),:), 3), 'b-', 'Linewidth', 2);
semilogx(K(1:K_max_index(3)), mean(bessel_enstrophy_adv_vel_alldirections_rmin_cutoff(3,1:K_max_index(3),:), 3), 'm-', 'Linewidth', 2);
semilogx(K(1:K_max_index(4)), mean(bessel_enstrophy_adv_vel_alldirections_rmin_cutoff(4,1:K_max_index(4),:), 3), 'r-', 'Linewidth', 2);

aomega = semilogx(K(1:K_max_index(1)), mean(bessel_enstrophy_adv_vort_alldirections_rmin_cutoff(1,1:K_max_index(1),:), 3), 'c:', 'Linewidth', 2);
semilogx(K(1:K_max_index(2)), mean(bessel_enstrophy_adv_vort_alldirections_rmin_cutoff(2,1:K_max_index(2),:), 3), 'b:', 'Linewidth', 2);
semilogx(K(1:K_max_index(3)), mean(bessel_enstrophy_adv_vort_alldirections_rmin_cutoff(3,1:K_max_index(3),:), 3), 'm:', 'Linewidth', 2);
semilogx(K(1:K_max_index(4)), mean(bessel_enstrophy_adv_vort_alldirections_rmin_cutoff(4,1:K_max_index(4),:), 3), 'r:', 'Linewidth', 2);

plot([K_max(1), K_max(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_max(2), K_max(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_max(3), K_max(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_max(4), K_max(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

ylabel('\Pi_K^{\omega}')

xlabel('Wavenumber K')
ylim([ymin, ymax]);
xlim([xmin, xmax]);

legend([Enstrophy_Flux aomega au], ...
    'Enstrophy Flux','SF_{A\omega}: Smallest r_{min}', 'SF_{Au}: Smallest r_{min}', ...
    'Location','NorthWest');
set(gca,'fontsize', size_of_font);

figure(21)

subplot(2,3,4)
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Enstrophy_Flux , 'k-', 'Linewidth', 5);
hold on

semilogx(K(K_min_index(4):end), mean(bessel_enstrophy_adv_vel_alldirections_rmax_cutoff(4,K_min_index(4):end,:), 3), 'c:', 'Linewidth', 1);
semilogx(K(K_min_index(3):end), mean(bessel_enstrophy_adv_vel_alldirections_rmax_cutoff(3,K_min_index(3):end,:), 3), 'b:', 'Linewidth', 1);
semilogx(K(K_min_index(2):end), mean(bessel_enstrophy_adv_vel_alldirections_rmax_cutoff(2,K_min_index(2):end,:), 3), 'm:', 'Linewidth', 1);
semilogx(K(K_min_index(1):end), mean(bessel_enstrophy_adv_vel_alldirections_rmax_cutoff(1,K_min_index(1):end,:), 3), 'r:', 'Linewidth', 1);
semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Enstrophy_Flux , 'k-', 'Linewidth', 5);

semilogx(K(K_min_index(1):end), mean(bessel_enstrophy_adv_vort_alldirections_rmax_cutoff(1,K_min_index(1):end,:), 3), 'c-', 'Linewidth', 2);
semilogx(K(K_min_index(2):end), mean(bessel_enstrophy_adv_vort_alldirections_rmax_cutoff(2,K_min_index(2):end,:), 3), 'b-', 'Linewidth', 2);
semilogx(K(K_min_index(3):end), mean(bessel_enstrophy_adv_vort_alldirections_rmax_cutoff(3,K_min_index(3):end,:), 3), 'm-', 'Linewidth', 2);
semilogx(K(K_min_index(4):end), mean(bessel_enstrophy_adv_vort_alldirections_rmax_cutoff(4,K_min_index(4):end,:), 3), 'r-', 'Linewidth', 2);

plot([K_min(1), K_min(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_min(2), K_min(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_min(3), K_min(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_min(4), K_min(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

ylabel('\Pi_K^{\omega}')

xlabel('Wavenumber K')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
set(gca,'fontsize', size_of_font);

%% Begin SQG Plotting

clear all

%% Start with SQG turbulence simulation

experiment_name = "StrongBeta"; %["SB","WB","SBd"];

quantiles = 0.75;

no_of_exps = size(experiment_name,2);

size_of_font = 14;
J_1_maximum = 2;
J_2_maximum = 3.1;
J_3_maximum = 4.3;
J_2_3_function_minimum = 1.5;
k_Halfbb_Cascade_End = 4e2;
r_Halfbb_Cascade_End = 2*pi/k_Halfbb_Cascade_End;

k_forcing_estimate = 1e1;
r_forcing_estimate = 2*pi/k_forcing_estimate;
data_path = "../analysis/processed_data/";

Diagonal = load(data_path+'Structure_Functions_SQG_'+experiment_name+'_Diagonal.mat');
Meridional = load(data_path+'Structure_Functions_SQG_'+experiment_name+'_Meridional.mat');
Zonal = load(data_path+'Structure_Functions_SQG_'+experiment_name+'_Zonal.mat');
OffDiagonal = load(data_path+'Structure_Functions_SQG_'+experiment_name+'_Off-Diagonal.mat');


%% Plot spectral fluxes

Spectral_Flux = load(data_path+'Spectral_Fluxes_SQG_'+experiment_name+'.mat');

Patch_k_Bounds = [Spectral_Flux.Wavenumber' fliplr(Spectral_Flux.Wavenumber')];
Patch_k_Bounds_Halfbb = [k_forcing_estimate, max(Spectral_Flux.Wavenumber), ...
    max(Spectral_Flux.Wavenumber), k_forcing_estimate];
%Patch_k_Bounds_Energy = [k_Energy_Cascade_End, k_forcing_estimate, ...
%    k_forcing_estimate, k_Energy_Cascade_End];

Spectral_Flux.bb_Flux_max = (quantile(Spectral_Flux.bb_Flux_snapshots ...
    ,quantiles,3))-Spectral_Flux.bb_Flux(1);
Spectral_Flux.bb_Flux_min = (quantile(Spectral_Flux.bb_Flux_snapshots ...
    ,1-quantiles,3))-Spectral_Flux.bb_Flux(1);
Spectral_Flux.KE_Flux_max = (quantile(Spectral_Flux.KE_Flux_snapshots ...
    ,quantiles,3));
Spectral_Flux.KE_Flux_min = (quantile(Spectral_Flux.KE_Flux_snapshots ...
    ,1-quantiles,3));

Spectral_Flux.bb_Flux = Spectral_Flux.bb_Flux/2;


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

n_min=[2; 4; 8; 16];
n_max = (1./n_min)*size(R_beta,1)/2;

n_min_D = ceil(n_min/sqrt(2));
n_max_D = floor(n_max/sqrt(2));


%% Calculate 2D turbulence Bessel function estimates of fluxes

R_spacing = R_beta(2)-R_beta(1);
R_spacing_Z = R_Intermediate(2)-R_Intermediate(1);

K = Spectral_Flux.Wavenumber;

for ii=1:size(K,1)
    % This loop iterates over wavenumber (K) values to find the flux that
    % would be estimated by various structure functions.
    % These estimates utilize Bessel function integrals over all
    % separation distances. Note that first run-through uses diagonal
    % ("Along_beta") separation vectors

    % Halfbb flux estimates. This example uses the advective 
    % bicity structure function 
   
    % First index represents the upper limit of 
    %  separation distance (R) used for the interval. Optimally should
    %  use the maximum R available, but for some purposes it is useful
    %  to know what flux you would estimate from an incomplete
    %  structure function dataset (only measuring up to a specific R)
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Along_beta.bbadv(:,1).*besselj(1,K(ii)*R_beta));

    % This is the Halfbb flux estimate
    bessel_Halfbb_adv_b_D(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_Halfbb_adv_b_rmin_cutoff_D(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_Halfbb_adv_b_rmax_cutoff_D(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end

    % Now repeat for energy fluxes using:

    % Advective velocity SF
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Along_beta.velveladv(:,1).*besselj(1,K(ii)*R_beta));
    bessel_energy_adv_vel_D(ii) = bessel_dummy_integral(end,ii);
    %bessel_energy_adv_vel_cutoff_D(ii) = bessel_dummy_integral(ii,ii);
    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_D(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_D(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
        r_min_D(k) = R_beta(n_min_D(k));
        r_max_D(k) = R_beta(n_max_D(k));
    end

    
    % Repeat all for zonal direction ("Intermediate_1")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.bbadv(:,1).*besselj(1,K(ii)*R_Intermediate));
    bessel_Halfbb_adv_b_Z(ii) = bessel_dummy_integral(end,ii);

     for k=1:size(n_min,1)
        bessel_Halfbb_adv_b_rmin_cutoff_Z(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_Halfbb_adv_b_rmax_cutoff_Z(k,ii) = bessel_dummy_integral(n_max(k),ii);
     end

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.velveladv(:,1).* ...
        besselj(1,K(ii)*R_Intermediate));
    bessel_energy_adv_vel_Z(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_Z(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_Z(k,ii) = bessel_dummy_integral(n_max(k),ii);
        r_min(k) = R_Intermediate(n_min(k));
        r_max(k) = R_Intermediate(n_max(k));
    end


    % Repeat all for meridional direction ("Intermediate_2")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.bbadv(:,1).*besselj(1,K(ii)*R_Intermediate));
    bessel_Halfbb_adv_b_M(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_Halfbb_adv_b_rmin_cutoff_M(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_Halfbb_adv_b_rmax_cutoff_M(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.velveladv(:,1).* ...
        besselj(1,K(ii)*R_Intermediate));
    bessel_energy_adv_vel_M(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_M(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_M(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end



    % Repeat all for off-diagonal direction ("Across_beta")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.bbadv(:,1).*besselj(1,K(ii)*R_beta));
    bessel_Halfbb_adv_b_OD(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_Halfbb_adv_b_rmin_cutoff_OD(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_Halfbb_adv_b_rmax_cutoff_OD(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end

    
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.velveladv(:,1).* ...
        besselj(1,K(ii)*R_beta));
    bessel_energy_adv_vel_OD(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_OD(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_OD(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end
    
end




%% Create simple plots for main paper text

%% Create Energy cascade plots

h7=figure(20)
set(h7,'Position',[10 10 1400 1400])
ymin = -2.5e-5;
ymax = 1e-5;
xmin = 0.7;
xmax = 300;

% Advective velocity SFs

K_max = J_1_maximum./r_min;
K_min = J_1_maximum./r_max;
for k=1:size(n_min,1)
    K_max_index(k) = size(K(K<K_max(k)),1);
    K_min_index(k) = size(K(K<K_min(k)),1); % Add 1 since first index used is the one after this limit
end

bessel_energy_adv_vel_alldirections = ...
    cat(3,bessel_energy_adv_vel_Z, bessel_energy_adv_vel_M, ...
    bessel_energy_adv_vel_D, bessel_energy_adv_vel_OD);

bessel_energy_adv_vel_alldirections_rmin_cutoff = ...
    cat(3,bessel_energy_adv_vel_rmin_cutoff_Z, bessel_energy_adv_vel_rmin_cutoff_M, ...
    bessel_energy_adv_vel_rmin_cutoff_D, bessel_energy_adv_vel_rmin_cutoff_OD);

bessel_energy_adv_vel_alldirections_rmax_cutoff = ...
    cat(3,bessel_energy_adv_vel_rmax_cutoff_Z, bessel_energy_adv_vel_rmax_cutoff_M, ...
    bessel_energy_adv_vel_rmax_cutoff_D, bessel_energy_adv_vel_rmax_cutoff_OD);

figure(20)
subplot(2,3,2)
KE_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux ,'k-', 'Linewidth', 5);
hold on

semilogx(K(1:K_max_index(1)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(1,1:K_max_index(1),:), 3), 'c-', 'Linewidth', 2);
semilogx(K(1:K_max_index(2)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(2,1:K_max_index(2),:), 3), 'b-', 'Linewidth', 2);
semilogx(K(1:K_max_index(3)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(3,1:K_max_index(3),:), 3), 'm-', 'Linewidth', 2);
semilogx(K(1:K_max_index(4)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(4,1:K_max_index(4),:), 3), 'r-', 'Linewidth', 2);

plot([K_max(1), K_max(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_max(2), K_max(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_max(3), K_max(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_max(4), K_max(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

%semilogx(J_1_maximum./R_beta,-0.5*Along_beta.velveladv(:,1),'r--', 'Linewidth', 1);

%plot([forcing_k, forcing_k], [ymin, ymax], 'k--', 'Linewidth', 2)
%plot([max_energy_cascade_k, max_energy_cascade_k], [ymin, ymax], 'k--', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

%semilogx(K, ...
 %   mean(bessel_energy_adv_vel_alldirections, 3),'r-', 'Linewidth', 2);
ylim([ymin, ymax]);
xlim([xmin, xmax]);
ylabel('\Pi_K^u')
%legend([KE_Flux Bessel_SFveladv_Flux], ...
%    'Actual Flux','Bessel Method: SF_{Au}', ...
%    'Location','NorthEast');
title("SQG Turbulence")
set(gca,'fontsize', size_of_font);

figure(21)
subplot(2,3,2)
KE_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux ,'k-', 'Linewidth', 5);
%Bessel_SFveladv_Flux = semilogx(K, ...
%    mean(bessel_energy_adv_vel_D, 3),'r-', 'Linewidth', 2);
hold on

semilogx(K(K_min_index(1):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(1,K_min_index(1):end,:), 3), 'c-', 'Linewidth', 2);
semilogx(K(K_min_index(2):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(2,K_min_index(2):end,:), 3), 'b-', 'Linewidth', 2);
semilogx(K(K_min_index(3):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(3,K_min_index(3):end,:), 3), 'm-', 'Linewidth', 2);
semilogx(K(K_min_index(4):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(4,K_min_index(4):end,:), 3), 'r-', 'Linewidth', 2);

plot([K_min(1), K_min(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_min(2), K_min(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_min(3), K_min(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_min(4), K_min(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

%semilogx(J_1_maximum./R_beta,-0.5*Along_beta.velveladv(:,1),'r--', 'Linewidth', 1);

%plot([forcing_k, forcing_k], [ymin, ymax], 'k--', 'Linewidth', 2)
%plot([max_energy_cascade_k, max_energy_cascade_k], [ymin, ymax], 'k--', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

%semilogx(K, ...
 %   mean(bessel_energy_adv_vel_alldirections, 3),'r-', 'Linewidth', 2);
ylim([ymin, ymax]);
xlim([xmin, xmax]);
ylabel('\Pi_K^u')
%legend([KE_Flux Bessel_SFveladv_Flux], ...
%    'Actual Flux','Bessel Method: SF_{Au}', ...
%    'Location','NorthEast');
title("SQG Turbulence")
set(gca,'fontsize', size_of_font);

%% Create Halfbb cascade plots

ymin = -8e-6;
ymax = 8e-6;

% Advective buoyancy SFs

bessel_Halfbb_adv_b_alldirections_rmin_cutoff = ...
    cat(3,bessel_Halfbb_adv_b_rmin_cutoff_Z, bessel_Halfbb_adv_b_rmin_cutoff_M, ...
    bessel_Halfbb_adv_b_rmin_cutoff_D, bessel_Halfbb_adv_b_rmin_cutoff_OD);

bessel_Halfbb_adv_b_alldirections_rmax_cutoff = ...
    cat(3,bessel_Halfbb_adv_b_rmax_cutoff_Z, bessel_Halfbb_adv_b_rmax_cutoff_M, ...
    bessel_Halfbb_adv_b_rmax_cutoff_D, bessel_Halfbb_adv_b_rmax_cutoff_OD);

figure(20)

subplot(2,3,5)
Halfbb_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.bb_Flux , 'k-', 'Linewidth', 5);
hold on

ab = semilogx(K(1:K_max_index(1)), mean(bessel_Halfbb_adv_b_alldirections_rmin_cutoff(1,1:K_max_index(1),:), 3), 'c-', 'Linewidth', 2);
semilogx(K(1:K_max_index(2)), mean(bessel_Halfbb_adv_b_alldirections_rmin_cutoff(2,1:K_max_index(2),:), 3), 'b-', 'Linewidth', 2);
semilogx(K(1:K_max_index(3)), mean(bessel_Halfbb_adv_b_alldirections_rmin_cutoff(3,1:K_max_index(3),:), 3), 'm-', 'Linewidth', 2);
semilogx(K(1:K_max_index(4)), mean(bessel_Halfbb_adv_b_alldirections_rmin_cutoff(4,1:K_max_index(4),:), 3), 'r-', 'Linewidth', 2);

plot([K_max(1), K_max(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_max(2), K_max(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_max(3), K_max(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_max(4), K_max(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)


xlabel('Wavenumber K')
ylabel('\Pi_K^b')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
legend([Halfbb_Flux ab], ...
    '0.5*Buoy. Var. Flux','SF_{Ab}: Smallest r_{min}', ...
    'Location','NorthWest');
set(gca,'fontsize', size_of_font);

figure(21)

subplot(2,3,5)
Halfbb_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.bb_Flux , 'k-', 'Linewidth', 5);
hold on

semilogx(K(K_min_index(1):end), mean(bessel_Halfbb_adv_b_alldirections_rmax_cutoff(1,K_min_index(1):end,:), 3), 'c-', 'Linewidth', 2);
semilogx(K(K_min_index(2):end), mean(bessel_Halfbb_adv_b_alldirections_rmax_cutoff(2,K_min_index(2):end,:), 3), 'b-', 'Linewidth', 2);
semilogx(K(K_min_index(3):end), mean(bessel_Halfbb_adv_b_alldirections_rmax_cutoff(3,K_min_index(3):end,:), 3), 'm-', 'Linewidth', 2);
semilogx(K(K_min_index(4):end), mean(bessel_Halfbb_adv_b_alldirections_rmax_cutoff(4,K_min_index(4):end,:), 3), 'r-', 'Linewidth', 2);

plot([K_min(1), K_min(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_min(2), K_min(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_min(3), K_min(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_min(4), K_min(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

xlabel('Wavenumber K')
ylabel('\Pi_K^b')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
set(gca,'fontsize', size_of_font);


%% Create QG panels

clear all

size_of_font = 14;
J_1_maximum = 2;
J_2_maximum = 3.1;
J_3_maximum = 4.3;
J_2_3_function_minimum = 1.5;

%% Start with QG turbulence simulation

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

n_min=[2; 4; 8; 16];
n_max = (1./n_min)*size(R_beta,1)/2;

n_min_D = ceil(n_min/sqrt(2));
n_max_D = floor(n_max/sqrt(2));


%% Calculate 2D turbulence Bessel function estimates of fluxes

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

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vort_rmin_cutoff_D(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_enstrophy_adv_vort_rmax_cutoff_D(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end

    % Now repeat for advective velocity SF estimate of enstrophy flux
    bessel_dummy_integral_D(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing,Along_beta.velveladv(:,1).* (2*besselj(2,K(ii)*R_beta)./(K(ii)*R_beta) - besselj(1,K(ii)*R_beta) ));
    bessel_enstrophy_adv_vel_D(ii) = bessel_dummy_integral(end,ii);
    
    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vel_rmin_cutoff_D(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_enstrophy_adv_vel_rmax_cutoff_D(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end
    


    % Now repeat for energy fluxes using:

    % Advective velocity SF
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Along_beta.velveladv(:,1).*besselj(1,K(ii)*R_beta));
    bessel_energy_adv_vel_D(ii) = bessel_dummy_integral(end,ii);
    %bessel_energy_adv_vel_cutoff_D(ii) = bessel_dummy_integral(ii,ii);
    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_D(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_D(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
        r_min_D(k) = R_beta(n_min_D(k));
        r_max_D(k) = R_beta(n_max_D(k));
    end

    
    % Repeat all for zonal direction ("Intermediate_1")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.vortvortadv(:,1).*besselj(1,K(ii)*R_Intermediate));
    bessel_enstrophy_adv_vort_Z(ii) = bessel_dummy_integral(end,ii);

     for k=1:size(n_min,1)
        bessel_enstrophy_adv_vort_rmin_cutoff_Z(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_enstrophy_adv_vort_rmax_cutoff_Z(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end

    bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing_Z,Intermediate_1.velveladv(:,1).* ...
        (2*besselj(2,K(ii)*R_Intermediate)./(K(ii)*R_Intermediate) - besselj(1,K(ii)*R_Intermediate) ));
    bessel_enstrophy_adv_vel_Z(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vel_rmin_cutoff_Z(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_enstrophy_adv_vel_rmax_cutoff_Z(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_1.velveladv(:,1).* ...
        besselj(1,K(ii)*R_Intermediate));
    bessel_energy_adv_vel_Z(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_Z(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_Z(k,ii) = bessel_dummy_integral(n_max(k),ii);
        r_min(k) = R_Intermediate(n_min(k));
        r_max(k) = R_Intermediate(n_max(k));
    end


    % Repeat all for meridional direction ("Intermediate_2")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.vortvortadv(:,1).*besselj(1,K(ii)*R_Intermediate));
    bessel_enstrophy_adv_vort_M(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vort_rmin_cutoff_M(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_enstrophy_adv_vort_rmax_cutoff_M(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end

    bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing_Z,Intermediate_2.velveladv(:,1).* ...
        (2*besselj(2,K(ii)*R_Intermediate)./(K(ii)*R_Intermediate) - besselj(1,K(ii)*R_Intermediate) ));
    bessel_enstrophy_adv_vel_M(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vel_rmin_cutoff_M(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_enstrophy_adv_vel_rmax_cutoff_M(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end
    
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing_Z,Intermediate_2.velveladv(:,1).* ...
        besselj(1,K(ii)*R_Intermediate));
    bessel_energy_adv_vel_M(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_M(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_M(k,ii) = bessel_dummy_integral(n_max(k),ii);
    end



    % Repeat all for off-diagonal direction ("Across_beta")

    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.vortvortadv(:,1).*besselj(1,K(ii)*R_beta));
    bessel_enstrophy_adv_vort_OD(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vort_rmin_cutoff_OD(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_enstrophy_adv_vort_rmax_cutoff_OD(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end

    bessel_dummy_integral(:,ii) = 0.5*(K(ii)^3)*cumtrapz(R_spacing,Across_beta.velveladv(:,1).* ...
        (2*besselj(2,K(ii)*R_beta)./(K(ii)*R_beta) - besselj(1,K(ii)*R_beta) ));
    bessel_enstrophy_adv_vel_OD(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_enstrophy_adv_vel_rmin_cutoff_OD(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_enstrophy_adv_vel_rmax_cutoff_OD(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end
    
    bessel_dummy_integral(:,ii) = -0.5*K(ii)*cumtrapz(R_spacing,Across_beta.velveladv(:,1).* ...
        besselj(1,K(ii)*R_beta));
    bessel_energy_adv_vel_OD(ii) = bessel_dummy_integral(end,ii);

    for k=1:size(n_min,1)
        bessel_energy_adv_vel_rmin_cutoff_OD(k,ii) = bessel_dummy_integral(end,ii) - bessel_dummy_integral(n_min_D(k),ii);
        bessel_energy_adv_vel_rmax_cutoff_OD(k,ii) = bessel_dummy_integral(n_max_D(k),ii);
    end
    
end




%% Create simple plots for main paper text

%% Create Energy cascade plots

h7=figure(20)
set(h7,'Position',[10 10 1400 1400])
ymin = -600;
ymax = 200;
xmin = 0.7;
xmax = 250;

% Advective velocity SFs

K_max = J_1_maximum./r_min;
K_min = J_1_maximum./r_max;
for k=1:size(n_min,1)
    K_max_index(k) = size(K(K<K_max(k)),1);
    K_min_index(k) = size(K(K<K_min(k)),1); % Add 1 since first index used is the one after this limit
end

bessel_energy_adv_vel_alldirections = ...
    cat(3,bessel_energy_adv_vel_Z, bessel_energy_adv_vel_M, ...
    bessel_energy_adv_vel_D, bessel_energy_adv_vel_OD);

bessel_energy_adv_vel_alldirections_rmin_cutoff = ...
    cat(3,bessel_energy_adv_vel_rmin_cutoff_Z, bessel_energy_adv_vel_rmin_cutoff_M, ...
    bessel_energy_adv_vel_rmin_cutoff_D, bessel_energy_adv_vel_rmin_cutoff_OD);

bessel_energy_adv_vel_alldirections_rmax_cutoff = ...
    cat(3,bessel_energy_adv_vel_rmax_cutoff_Z, bessel_energy_adv_vel_rmax_cutoff_M, ...
    bessel_energy_adv_vel_rmax_cutoff_D, bessel_energy_adv_vel_rmax_cutoff_OD);

figure(20)
subplot(2,3,3)
KE_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux(:,1) ,'k-', 'Linewidth', 5);
hold on

semilogx(K(1:K_max_index(1)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(1,1:K_max_index(1),:), 3), 'c-', 'Linewidth', 2);
semilogx(K(1:K_max_index(2)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(2,1:K_max_index(2),:), 3), 'b-', 'Linewidth', 2);
semilogx(K(1:K_max_index(3)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(3,1:K_max_index(3),:), 3), 'm-', 'Linewidth', 2);
semilogx(K(1:K_max_index(4)), mean(bessel_energy_adv_vel_alldirections_rmin_cutoff(4,1:K_max_index(4),:), 3), 'r-', 'Linewidth', 2);

plot([K_max(1), K_max(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_max(2), K_max(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_max(3), K_max(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_max(4), K_max(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

%semilogx(J_1_maximum./R_beta,-0.5*Along_beta.velveladv(:,1),'r--', 'Linewidth', 1);

%plot([forcing_k, forcing_k], [ymin, ymax], 'k--', 'Linewidth', 2)
%plot([max_energy_cascade_k, max_energy_cascade_k], [ymin, ymax], 'k--', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

%semilogx(K, ...
 %   mean(bessel_energy_adv_vel_alldirections, 3),'r-', 'Linewidth', 2);
xlabel('Wavenumber K')
ylabel('\Pi_K^u')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
%legend([KE_Flux Bessel_SFveladv_Flux], ...
%    'Actual Flux','Bessel Method: SF_{Au}', ...
%    'Location','NorthEast');
title("QG Turbulence")
set(gca,'fontsize', size_of_font);

figure(21)
subplot(2,3,3)
KE_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux(:,1) ,'k-', 'Linewidth', 5);
%Bessel_SFveladv_Flux = semilogx(K, ...
%    mean(bessel_energy_adv_vel_D, 3),'r-', 'Linewidth', 2);
hold on

semilogx(K(K_min_index(1):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(1,K_min_index(1):end,:), 3), 'c-', 'Linewidth', 2);
semilogx(K(K_min_index(2):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(2,K_min_index(2):end,:), 3), 'b-', 'Linewidth', 2);
semilogx(K(K_min_index(3):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(3,K_min_index(3):end,:), 3), 'm-', 'Linewidth', 2);
semilogx(K(K_min_index(4):end), mean(bessel_energy_adv_vel_alldirections_rmax_cutoff(4,K_min_index(4):end,:), 3), 'r-', 'Linewidth', 2);

plot([K_min(1), K_min(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_min(2), K_min(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_min(3), K_min(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_min(4), K_min(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

%semilogx(J_1_maximum./R_beta,-0.5*Along_beta.velveladv(:,1),'r--', 'Linewidth', 1);

%plot([forcing_k, forcing_k], [ymin, ymax], 'k--', 'Linewidth', 2)
%plot([max_energy_cascade_k, max_energy_cascade_k], [ymin, ymax], 'k--', 'Linewidth', 2)
plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

%semilogx(K, ...
 %   mean(bessel_energy_adv_vel_alldirections, 3),'r-', 'Linewidth', 2);
xlabel('Wavenumber K')
ylabel('\Pi_K^u')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
%legend([KE_Flux Bessel_SFveladv_Flux], ...
%    'Actual Flux','Bessel Method: SF_{Au}', ...
%    'Location','NorthEast');
title("QG Turbulence")
set(gca,'fontsize', size_of_font);

%% Create enstrophy cascade plots

ymin = -2e5;
ymax = 5e5;

% Advective vorticity SFs

bessel_enstrophy_adv_vel_alldirections = ...
    cat(3,bessel_enstrophy_adv_vel_Z, bessel_enstrophy_adv_vel_M, ...
    bessel_enstrophy_adv_vel_D, bessel_enstrophy_adv_vel_OD);

bessel_enstrophy_adv_vel_alldirections_rmin_cutoff = ...
    cat(3,bessel_enstrophy_adv_vel_rmin_cutoff_Z, bessel_enstrophy_adv_vel_rmin_cutoff_M, ...
    bessel_enstrophy_adv_vel_rmin_cutoff_D, bessel_enstrophy_adv_vel_rmin_cutoff_OD);

bessel_enstrophy_adv_vel_alldirections_rmax_cutoff = ...
    cat(3,bessel_enstrophy_adv_vel_rmax_cutoff_Z, bessel_enstrophy_adv_vel_rmax_cutoff_M, ...
    bessel_enstrophy_adv_vel_rmax_cutoff_D, bessel_enstrophy_adv_vel_rmax_cutoff_OD);

bessel_enstrophy_adv_vort_alldirections_rmin_cutoff = ...
    cat(3,bessel_enstrophy_adv_vort_rmin_cutoff_Z, bessel_enstrophy_adv_vort_rmin_cutoff_M, ...
    bessel_enstrophy_adv_vort_rmin_cutoff_D, bessel_enstrophy_adv_vort_rmin_cutoff_OD);

bessel_enstrophy_adv_vort_alldirections_rmax_cutoff = ...
    cat(3,bessel_enstrophy_adv_vort_rmax_cutoff_Z, bessel_enstrophy_adv_vort_rmax_cutoff_M, ...
    bessel_enstrophy_adv_vort_rmax_cutoff_D, bessel_enstrophy_adv_vort_rmax_cutoff_OD);

figure(20)
subplot(2,3,6)
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.vortvort_Flux(:,1) , 'k-', 'Linewidth', 5);
hold on

au = semilogx(K(1:K_max_index(1)), mean(bessel_enstrophy_adv_vel_alldirections_rmin_cutoff(1,1:K_max_index(1),:), 3), 'c:', 'Linewidth', 2);
semilogx(K(1:K_max_index(2)), mean(bessel_enstrophy_adv_vel_alldirections_rmin_cutoff(2,1:K_max_index(2),:), 3), 'b:', 'Linewidth', 2);
semilogx(K(1:K_max_index(3)), mean(bessel_enstrophy_adv_vel_alldirections_rmin_cutoff(3,1:K_max_index(3),:), 3), 'm:', 'Linewidth', 2);
semilogx(K(1:K_max_index(4)), mean(bessel_enstrophy_adv_vel_alldirections_rmin_cutoff(4,1:K_max_index(4),:), 3), 'r:', 'Linewidth', 2);

aomega = semilogx(K(1:K_max_index(1)), mean(bessel_enstrophy_adv_vort_alldirections_rmin_cutoff(1,1:K_max_index(1),:), 3), 'c-', 'Linewidth', 2);
semilogx(K(1:K_max_index(2)), mean(bessel_enstrophy_adv_vort_alldirections_rmin_cutoff(2,1:K_max_index(2),:), 3), 'b-', 'Linewidth', 2);
semilogx(K(1:K_max_index(3)), mean(bessel_enstrophy_adv_vort_alldirections_rmin_cutoff(3,1:K_max_index(3),:), 3), 'm-', 'Linewidth', 2);
semilogx(K(1:K_max_index(4)), mean(bessel_enstrophy_adv_vort_alldirections_rmin_cutoff(4,1:K_max_index(4),:), 3), 'r-', 'Linewidth', 2);

plot([K_max(1), K_max(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_max(2), K_max(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_max(3), K_max(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_max(4), K_max(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)


xlabel('Wavenumber K')
ylabel('\Pi_K^{\omega}')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
legend([Enstrophy_Flux aomega au], ...
    'Enstrophy Flux','SF_{A\omega}: Smallest r_{min}', 'SF_{Au}: Smallest r_{min}', ...
    'Location','NorthWest');
set(gca,'fontsize', size_of_font);


figure(21)
subplot(2,3,6)
Enstrophy_Flux = semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.vortvort_Flux(:,1) , 'k-', 'Linewidth', 5);
hold on

semilogx(K(K_min_index(1):end), mean(bessel_enstrophy_adv_vel_alldirections_rmax_cutoff(1,K_min_index(1):end,:), 3), 'c:', 'Linewidth', 1);
semilogx(K(K_min_index(2):end), mean(bessel_enstrophy_adv_vel_alldirections_rmax_cutoff(2,K_min_index(2):end,:), 3), 'b:', 'Linewidth', 1);
semilogx(K(K_min_index(3):end), mean(bessel_enstrophy_adv_vel_alldirections_rmax_cutoff(3,K_min_index(3):end,:), 3), 'm:', 'Linewidth', 1);
semilogx(K(K_min_index(4):end), mean(bessel_enstrophy_adv_vel_alldirections_rmax_cutoff(4,K_min_index(4):end,:), 3), 'r:', 'Linewidth', 1);

semilogx(K(K_min_index(1):end), mean(bessel_enstrophy_adv_vort_alldirections_rmax_cutoff(1,K_min_index(1):end,:), 3), 'c-', 'Linewidth', 2);
semilogx(K(K_min_index(2):end), mean(bessel_enstrophy_adv_vort_alldirections_rmax_cutoff(2,K_min_index(2):end,:), 3), 'b-', 'Linewidth', 2);
semilogx(K(K_min_index(3):end), mean(bessel_enstrophy_adv_vort_alldirections_rmax_cutoff(3,K_min_index(3):end,:), 3), 'm-', 'Linewidth', 2);
semilogx(K(K_min_index(4):end), mean(bessel_enstrophy_adv_vort_alldirections_rmax_cutoff(4,K_min_index(4):end,:), 3), 'r-', 'Linewidth', 2);

plot([K_min(1), K_min(1)], [ymin, ymax], 'c--', 'Linewidth', 2)
plot([K_min(2), K_min(2)], [ymin, ymax], 'b--', 'Linewidth', 2)
plot([K_min(3), K_min(3)], [ymin, ymax], 'm--', 'Linewidth', 2)
plot([K_min(4), K_min(4)], [ymin, ymax], 'r--', 'Linewidth', 2)

plot([xmin, xmax], [0, 0], 'k-', 'Linewidth', 1)

xlabel('Wavenumber K')
ylabel('\Pi_K^{\omega}')
ylim([ymin, ymax]);
xlim([xmin, xmax]);
set(gca,'fontsize', size_of_font);
