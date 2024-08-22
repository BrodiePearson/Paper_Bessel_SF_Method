%% Code to diagnose spectral fluxes maps for each layer from pyqg data and
%  their isotropic averages by using spectral transforms consistent with
%  pyqg.

clear all

% Define the beta values for the set of 4 experiments
% Strong meridional beta (SB)
% Strong diagonal beta (SBd)
% Weak meridional beta (WB)
% Weak diagonal beta (WBd)
% No beta (NB)
experiment = ["SB","SBd","WB","WBd","NB"];
betay = ["10.0", "7.071067811865475", "1.0", "0.7071067811865475", "0.0"];
betax = ["0.0", "7.071067811865475", "0.0", "0.7071067811865475", "0.0"];

for exp_no = 1:size(experiment,2)

filedir = "../simulations/Output/Data/Equilibrated/";
files = dir(filedir + "Anisotropic2D_n_2048_drag_0.04_order_-2_visc"+ ...
    "_1.0e-21_order_8_kf_100.0_F_1.0e-5_betay_"+betay(exp_no)+"_"+ ...
    "betax_"+betax(exp_no)+"_*.mat");

% Sort the snapshots from earliest to latest
temp1 = struct2table(files); % convert the struct array to a table
temp2 = sortrows(temp1, 'datenum'); % sort the table by 'DOB'
files = table2struct(temp2);

clear temp1 temp2

file_no = length(files);

%% Load a data file to set size of structure function arrays before entering
%  parallel loop

S1 = load(filedir + files(1).name);
q = S1.zeta;
x = [S1.dx/2:S1.dx:S1.Lx-S1.dx/2];
y = [S1.dy/2:S1.dy:S1.Lx-S1.dy/2];
u=S1.u;
v=S1.v;
% Tile x and y positions onto grid
x=repmat(x',1,size(x,2));
y=repmat(y,size(y,2),1);

domain_width = 2*pi; % Domain width and length (square)
grid_width = x(2,1)-x(1,1); % Zonal and meridional grid-spacing

dX = grid_width;    % Spatial sampling interval
N = size(x,1);   % Number of spatial points in each direction (isotropic)
k_int = 1/dX;    % Smallest wavenumber [1/s] NOT [rad/s]
k = 2*pi*(-N/2:N/2-1)*(k_int/N); % x-wavenumber
l = k; % y-wavenumber (isotropic)

k_max_mat = max(max(-k(1:end/2+1)),max(l));
dk = 2*pi/domain_width;
dl=dk;
dkr_mat = sqrt(dk.^2+dl.^2);
kr = dkr_mat/2:dkr_mat:k_max_mat+dkr_mat;
kr=kr';

kr_size = size(kr,1);

clearvars -except kr_size file_no filedir files exp_no experiment betay betax


%% Calculate structure functions for each snapshot using parallel loops

parfor nn=1:file_no
    
    tmp_Enstrophy_Flux=zeros(kr_size,1);
    tmp_Energy_Flux=zeros(kr_size,1);
    
    S = load(filedir + files(nn).name)
    q = S.zeta;
    x = [S.dx/2:S.dx:S.Lx-S.dx/2];
    y = [S.dy/2:S.dy:S.Lx-S.dy/2];
    u = S.u;
    v = S.v;
    Work_Rate_ENS = S.Work_Rate_ENS;
    Drag_Rate_ENS = S.Drag_Rate_ENS;
    Diss_Rate_ENS = S.Diss_Rate_ENS;
    Work_Rate_KE = S.Work_Rate_KE;
    Drag_Rate_KE = S.Drag_Rate_KE;
    Diss_Rate_KE = S.Diss_Rate_KE;
    Energy = S.KE;
    % Tile x and y positions onto grid
    x=repmat(x',1,size(x,2));
    y=repmat(y,size(y,2),1);
    
    domain_width = 2*pi; % Domain width and length (square)
    grid_width = x(2,1)-x(1,1); % Zonal and meridional grid-spacing
    
    %% Calculate gradients and advective terms for new structure functions
    
    % Calculate gradients in velocity (for advection and vorticity)
    % Note that for a variable f(i,j), i is the y-component and j is the
    % x-component
    x_sep=abs(x(2,1)-x(1,1));
    y_sep=abs(y(1,2)-y(1,1));
    dudx=(circshift(u,-1,1)-circshift(u,1,1))./(2*x_sep);
    dudy=(circshift(u,-1,2)-circshift(u,1,2))./(2*y_sep);
    dvdx=(circshift(v,-1,1)-circshift(v,1,1))./(2*x_sep);
    dvdy=(circshift(v,-1,2)-circshift(v,1,2))./(2*y_sep);
    
    q = dudy-dvdx;
    
    % Calculate true (q) vorticity advection
    dqdx=(circshift(q,-1,1)-circshift(q,1,1))./(2*x_sep);
    dqdy=(circshift(q,-1,2)-circshift(q,1,2))./(2*y_sep);
    
    %% Diagnose Spectral Flux of Enstrophy and KE
    
    dX = grid_width;    % Spatial sampling interval
    N = size(x,1);   % Number of spatial points in each direction (isotropic)
    k_int = 1/dX;    % Smallest wavenumber [1/s] NOT [rad/s]
    k = 2*pi*(-N/2:N/2-1)*(k_int/N); % x-wavenumber
    l = k; % y-wavenumber (isotropic)
    [k_mat, l_mat] = meshgrid(k,l); % 2D-field of k and l for matlab transforms
    
    % Calculate gradients exactly via spectral calculation
    dqdx = real(ifft2(ifftshift(1i*k_mat'.*fftshift(fft2(q)))));
    dqdy = real(ifft2(ifftshift(1i*l_mat'.*fftshift(fft2(q)))));
    dudx = real(ifft2(ifftshift(1i*k_mat'.*fftshift(fft2(u)))));
    dudy = real(ifft2(ifftshift(1i*l_mat'.*fftshift(fft2(u)))));
    dvdx = real(ifft2(ifftshift(1i*k_mat'.*fftshift(fft2(v)))));
    dvdy = real(ifft2(ifftshift(1i*l_mat'.*fftshift(fft2(v)))));
    
    % Diagnose Fourier transformed vorticity advection by multiplying velocity
    % and vorticity in real space, Fourier transforming, multiplying by
    % the grid spacing (to convert fft's discrete sum to a discrete integral
    % over space) and dividing by 2\pi to convert to a radial wavenumber
    % transform
    J_q=fft2(u.*dqdx + v.*dqdy)*dX^2/(2*pi);
    J_u=fft2(u.*dudx + v.*dudy)*dX^2/(2*pi);
    J_v=fft2(u.*dvdx + v.*dvdy)*dX^2/(2*pi);
    
    q_f = fft2(q)*dX^2/(2*pi);
    u_f = fft2(u)*dX^2/(2*pi);
    v_f = fft2(v)*dX^2/(2*pi);
    
    % Derive the spectral flux divergences from advection and streamfunction
    ens_flux_div = fftshift(real(conj(q_f).*J_q));
    KE_flux_div = fftshift(real(conj(u_f).*J_u + conj(v_f).*J_v));
    [k_mat, l_mat] = meshgrid(k,l);
    
    k_max_mat = max(max(-k(1:end/2+1)),max(l));
    dk = 2*pi/domain_width;
    dl=dk;
    dkr_mat = sqrt(dk.^2+dl.^2);
    tmp_kr = dkr_mat/2:dkr_mat:k_max_mat+dkr_mat;
    tmp_kr=tmp_kr';
    
    % Calculate Fluxes by summing the flux divergences over
    for kk=1:size(tmp_kr,1)
        tmp_Enstrophy_Flux(kk) = -sum(ens_flux_div(k_mat.^2 + ...
            l_mat.^2 >= tmp_kr(kk)^2)*dk*dl)/((2*pi)^2);
        tmp_Energy_Flux(kk) = -sum(KE_flux_div(k_mat.^2 + ...
            l_mat.^2 >= tmp_kr(kk)^2)*dk*dl)/((2*pi)^2);
    end
    
    %% Calculate Spectral Fluxes from data
    
    Enstrophy_Flux(:,nn) = tmp_Enstrophy_Flux
    Energy_Flux(:,nn) = tmp_Energy_Flux
    kr(:,nn) = tmp_kr;
    
    Enstrophy_Work_Rate(nn,1) = Work_Rate_ENS;
    Enstrophy_Drag_Rate(nn,1) = Drag_Rate_ENS;
    Enstrophy_Diss_Rate(nn,1) = Diss_Rate_ENS;
    Energy_Work_Rate(nn,1) = Work_Rate_KE;
    Energy_Drag_Rate(nn,1) = Drag_Rate_KE;
    Energy_Diss_Rate(nn,1) = Diss_Rate_KE;
    time(nn)=files(nn).datenum;
    
end

%% Save structure functions averaged across snapshots to a structure

Spectral_Flux.Energy_Flux = nanmean(Energy_Flux,2);
Spectral_Flux.Enstrophy_Flux = nanmean(Enstrophy_Flux,2);

Spectral_Flux.Energy_Flux_snapshots = Energy_Flux;
Spectral_Flux.Enstrophy_Flux_snapshots = Enstrophy_Flux;

Spectral_Flux.Wavenumber = kr(:,1);

%% Save enstrophy and energy diagnostics to structure

Spectral_Flux.Enstrophy_Work_Rate = Enstrophy_Work_Rate;
Spectral_Flux.Enstrophy_Drag_Rate = Enstrophy_Drag_Rate;
Spectral_Flux.Enstrophy_Diss_Rate = Enstrophy_Diss_Rate;

Spectral_Flux.Energy_Work_Rate = Energy_Work_Rate;
Spectral_Flux.Energy_Drag_Rate = Energy_Drag_Rate;
Spectral_Flux.Energy_Diss_Rate = Energy_Diss_Rate;
Spectral_Flux.Time_of_snapshots = time;

%% Save Structure Functions and their errors into a .mat file

save("Spectral_Fluxes_"+experiment(exp_no)+'discretization_check','-struct','Spectral_Flux','-v7.3');


%% Plot spectral fluxes in spectral space

figure
subplot(1,2,1)
semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Enstrophy_Flux,'k-')
hold on
semilogx(Spectral_Flux.Wavenumber, ...
    Spectral_Flux.Enstrophy_Flux+nanstd(Spectral_Flux.Enstrophy_Flux_snapshots,0,2),'k:')
semilogx(Spectral_Flux.Wavenumber, ...
    Spectral_Flux.Enstrophy_Flux-nanstd(Spectral_Flux.Enstrophy_Flux_snapshots,0,2),'k:')
semilogx([100 max(Spectral_Flux.Wavenumber)], ...
    Spectral_Flux.Enstrophy_Diss_Rate*[1 1], 'r--')
semilogx([min(Spectral_Flux.Wavenumber) 100], ...
    -Spectral_Flux.Enstrophy_Drag_Rate*[1 1], 'r--')
semilogx([100 100], ...
    [-Spectral_Flux.Enstrophy_Drag_Rate Spectral_Flux.Enstrophy_Diss_Rate], 'r--')
title('Enstrophy Flux: experiment'+experiment(exp_no))
subplot(1,2,2)
semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.Energy_Flux,'k-')
hold on
semilogx(Spectral_Flux.Wavenumber, ...
    Spectral_Flux.Energy_Flux+nanstd(Spectral_Flux.Energy_Flux_snapshots,0,2),'k:')
semilogx(Spectral_Flux.Wavenumber, ...
    Spectral_Flux.Energy_Flux-nanstd(Spectral_Flux.Energy_Flux_snapshots,0,2),'k:')
semilogx([100 max(Spectral_Flux.Wavenumber)], ...
    Spectral_Flux.Energy_Diss_Rate*[1 1], 'r--')
semilogx([min(Spectral_Flux.Wavenumber) 100], ...
    -Spectral_Flux.Energy_Drag_Rate*[1 1], 'r--')
semilogx([100 100], ...
    [-Spectral_Flux.Energy_Drag_Rate Spectral_Flux.Energy_Diss_Rate], 'r--')
title('KE Flux')

end

