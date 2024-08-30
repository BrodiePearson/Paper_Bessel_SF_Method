%% Code to diagnose Multi-layer QG spectral fluxes from GeophysicalFlows.jl

clear all

nx = "256";
hyperviscosity = "1.0e-17";
hyperviscosity_order = "8";
hypoviscosity = "0.005";
hypoviscosity_order = "-2";
forcing_wavenumber = "7.0";
forcing_rate = "1.0e-5";
endtime = "30000.0";
% beta = "0.5"; % Switch beta off for some simulations where it is zero
% Let the model spin up for two eddy turnover times (one eddy turnover time ~500)
% So remove files that are before time = 1000
start_file = 17; %Time = 40 or 2/mu which is when the simulation equilibrates

experiment_name = "multilayerqg_2layer_" + nx ;

filedir = "../simulations/QG_simulations/data_2layer/";

files = dir(filedir + experiment_name +"_*.mat");

% Sort the snapshots from earliest to latest
temp1 = struct2table(files); % convert the struct array to a table
temp2 = sortrows(temp1, 'datenum'); % sort the table by 'DOB'
files = table2struct(temp2);

clear temp1 temp2

files(1:start_file) = [];
file_no = length(files);

%% Load a data file to set size of structure function arrays before entering
%  parallel loop

S1 = load(filedir + files(1).name);
q = S1.q; % vorticity
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

clearvars -except kr_size file_no filedir files experiment_name


%% Calculate structure functions for each snapshot using parallel loops

for layer=1:2
    for nn=1:file_no

        tmp_qq_Flux=zeros(kr_size,1);
        tmp_KE_Flux=zeros(kr_size,1);

        S = load(filedir + files(nn).name);
        q = S.q(:,:,layer);
        x = [S.dx/2:S.dx:S.Lx-S.dx/2];
        y = [S.dy/2:S.dy:S.Lx-S.dy/2];
        u = S.u(:,:,layer);
        v = S.v(:,:,layer);
        % buoyancy_work = S.buoyancy_work;
        % buoyancy_dissipation_hypoviscosity = S.buoyancy_dissipation_hypoviscosity;
        % buoyancy_dissipation_hyperviscosity = S.buoyancy_dissipation_hyperviscosity;
        KE = S.KE(layer);
        % Tile x and y positions onto grid
        x=repmat(x',1,size(x,2));
        y=repmat(y,size(y,2),1);

        domain_width = 2*pi; % Domain width and length (square)
        grid_width = x(2,1)-x(1,1); % Zonal and meridional grid-spacing

        %% Diagnose spatial/spectral domain variables

        dX = grid_width;    % Spatial sampling interval
        N = size(x,1);   % Number of spatial points in each direction (isotropic)
        k_int = 1/dX;    % Smallest wavenumber [1/s] NOT [rad/s]
        k = 2*pi*(-N/2:N/2-1)*(k_int/N); % x-wavenumber
        l = k; % y-wavenumber (isotropic)
        [k_mat, l_mat] = meshgrid(k,l); % 2D-field of k and l for matlab transforms

        %% Calculate gradients and advective terms for new structure functions

        % Calculate gradients exactly via spectral calculation
        dqdx = real(ifft2(ifftshift(1i*k_mat'.*fftshift(fft2(q)))));
        dqdy = real(ifft2(ifftshift(1i*l_mat'.*fftshift(fft2(q)))));
        dudx = real(ifft2(ifftshift(1i*k_mat'.*fftshift(fft2(u)))));
        dudy = real(ifft2(ifftshift(1i*l_mat'.*fftshift(fft2(u)))));
        dvdx = real(ifft2(ifftshift(1i*k_mat'.*fftshift(fft2(v)))));
        dvdy = real(ifft2(ifftshift(1i*l_mat'.*fftshift(fft2(v)))));

        % Calculate relative vorticity and its gradients via spectral calc.
        vort = dvdx-dudy;
        dvortdx = real(ifft2(ifftshift(1i*k_mat'.*fftshift(fft2(vort)))));
        dvortdy = real(ifft2(ifftshift(1i*l_mat'.*fftshift(fft2(vort)))));


        % Calculate gradients using finite differencing (_fd)
        % Note that for a variable f(i,j), i is the y-component and j is the
        % x-component
        x_sep=abs(x(2,1)-x(1,1));
        y_sep=abs(y(1,2)-y(1,1));

        dudx_fd=(circshift(u,-1,1)-circshift(u,1,1))./(2*x_sep);
        dudy_fd=(circshift(u,-1,2)-circshift(u,1,2))./(2*y_sep);
        dvdx_fd=(circshift(v,-1,1)-circshift(v,1,1))./(2*x_sep);
        dvdy_fd=(circshift(v,-1,2)-circshift(v,1,2))./(2*y_sep);
        dqdx_fd=(circshift(q,-1,1)-circshift(q,1,1))./(2*x_sep);
        dqdy_fd=(circshift(q,-1,2)-circshift(q,1,2))./(2*y_sep);

        % Calculate relative vorticity and its gradients via fd
        vort_fd = dvdx_fd-dudy_fd;
        dvortdx_fd = (circshift(vort_fd,-1,1)-circshift(vort_fd,1,1))./(2*x_sep);
        dvortdy_fd = (circshift(vort_fd,-1,2)-circshift(vort_fd,1,2))./(2*y_sep);


        %% Diagnose Spectral Flux of PE and KE

        % Diagnose Fourier transformed buoyancy advection by multiplying velocity
        % and vorticity in real space, Fourier transforming, multiplying by
        % the grid spacing (to convert fft's discrete sum to a discrete integral
        % over space) and dividing by 2\pi to convert to a radial wavenumber
        % transform
        J_q=fft2(u.*dqdx + v.*dqdy)*dX^2/(2*pi);
        J_u=fft2(u.*dudx + v.*dudy)*dX^2/(2*pi);
        J_v=fft2(u.*dvdx + v.*dvdy)*dX^2/(2*pi);
        J_vort=fft2(u.*dvortdx + v.*dvortdy)*dX^2/(2*pi);

        J_q_fd=fft2(u.*dqdx_fd + v.*dqdy_fd)*dX^2/(2*pi);
        J_u_fd=fft2(u.*dudx_fd + v.*dudy_fd)*dX^2/(2*pi);
        J_v_fd=fft2(u.*dvdx_fd + v.*dvdy_fd)*dX^2/(2*pi);
        J_vort_fd=fft2(u.*dvortdx_fd + v.*dvortdy_fd)*dX^2/(2*pi);

        q_f = fft2(q)*dX^2/(2*pi);
        u_f = fft2(u)*dX^2/(2*pi);
        v_f = fft2(v)*dX^2/(2*pi);
        vort_f = fft2(vort)*dX^2/(2*pi);

        % Derive the spectral flux divergences from advection and streamfunction
        qq_flux_div = fftshift(real(conj(q_f).*J_q));
        KE_flux_div = fftshift(real(conj(u_f).*J_u + conj(v_f).*J_v));
        vortvort_flux_div = fftshift(real(conj(vort_f).*J_vort));

        qq_flux_div_fd = fftshift(real(conj(q_f).*J_q_fd));
        KE_flux_div_fd = fftshift(real(conj(u_f).*J_u_fd + conj(v_f).*J_v_fd));
        vortvort_flux_div_fd = fftshift(real(conj(vort_f).*J_vort_fd));

        qq_spectra = fftshift(real(conj(q_f).*q_f));
        vortvort_spectra = fftshift(real(conj(vort_f).*vort_f));

        k_max_mat = max(max(-k(1:end/2+1)),max(l));
        dk = 2*pi/domain_width;
        dl=dk;
        dkr_mat = sqrt(dk.^2+dl.^2);
        tmp_kr = dkr_mat/2:dkr_mat:k_max_mat+dkr_mat;
        tmp_kr=tmp_kr';

        % Calculate Fluxes of buoyancy variance and kinetic energy
        % by summing the flux divergences over all wavenumbers
        for kk=1:size(tmp_kr,1)
            tmp_qq_Flux(kk) = -sum(qq_flux_div(k_mat.^2 + ...
                l_mat.^2 >= tmp_kr(kk)^2)*dk*dl)/((2*pi)^2);
            tmp_KE_Flux(kk) = -sum(KE_flux_div(k_mat.^2 + ...
                l_mat.^2 >= tmp_kr(kk)^2)*dk*dl)/((2*pi)^2);
            tmp_vortvort_Flux(kk) = -sum(vortvort_flux_div(k_mat.^2 + ...
                l_mat.^2 >= tmp_kr(kk)^2)*dk*dl)/((2*pi)^2);

            tmp_qq_Flux_fd(kk) = -sum(qq_flux_div_fd(k_mat.^2 + ...
                l_mat.^2 >= tmp_kr(kk)^2)*dk*dl)/((2*pi)^2);
            tmp_KE_Flux_fd(kk) = -sum(KE_flux_div_fd(k_mat.^2 + ...
                l_mat.^2 >= tmp_kr(kk)^2)*dk*dl)/((2*pi)^2);
            tmp_vortvort_Flux_fd(kk) = -sum(vortvort_flux_div_fd(k_mat.^2 + ...
                l_mat.^2 >= tmp_kr(kk)^2)*dk*dl)/((2*pi)^2);

            %tmp_qq_spectrum = qq_spectraqq_spectra = fftshift(real(conj(b_f).*b_f)); = fftshift(real(conj(b_f).*b_f));
        end

        %% Calculate Spectral Fluxes from data

        % Shift the spectra so zero at lowest wavenumber/largest scale
        qq_Flux(:,layer,nn) = tmp_qq_Flux ;%- tmp_qq_Flux(1);
        KE_Flux(:,layer,nn) = tmp_KE_Flux ;%- tmp_KE_Flux(1);
        vortvort_Flux(:,layer,nn) = tmp_vortvort_Flux ;%- tmp_qq_Flux(1);

        qq_Flux_fd(:,layer,nn) = tmp_qq_Flux_fd ;%- tmp_qq_Flux(1);
        KE_Flux_fd(:,layer,nn) = tmp_KE_Flux_fd ;%- tmp_KE_Flux(1);
        vortvort_Flux_fd(:,layer,nn) = tmp_vortvort_Flux_fd ;%- tmp_qq_Flux(1);
        qq_Flux_fd_flip(:,layer,nn) = tmp_qq_Flux_fd - tmp_qq_Flux_fd(1);
        KE_Flux_fd_flip(:,layer,nn) = tmp_KE_Flux_fd - tmp_KE_Flux_fd(1);
        vortvort_Flux_fd_flip(:,layer,nn) = tmp_vortvort_Flux_fd - tmp_vortvort_Flux_fd(1);

        kr(:,layer,nn) = tmp_kr;
        kinetic_energy(layer,nn) = KE;
        %qq_spectrum(:,nn) = qq_spectra; %- tmp_qq_Flux(1);

        % qq_Work_Rate(nn,1) = buoyancy_work;
        % qq_Drag_Rate(nn,1) = buoyancy_dissipation_hypoviscosity;
        % qq_Diss_Rate(nn,1) = buoyancy_dissipation_hyperviscosity;
        time(layer,nn)=files(nn).datenum;

    end
end

%% Save structure functions averaged across snapshots to a structure

Spectral_Flux.KE_Flux = mean(KE_Flux,3, 'omitnan');
Spectral_Flux.qq_Flux = mean(qq_Flux,3, 'omitnan');
Spectral_Flux.vortvort_Flux = mean(vortvort_Flux,3, 'omitnan');

Spectral_Flux.KE_Flux_fd = mean(KE_Flux_fd,3, 'omitnan');
Spectral_Flux.qq_Flux_fd = mean(qq_Flux_fd,3, 'omitnan');
Spectral_Flux.vortvort_Flux_fd = mean(vortvort_Flux_fd,3, 'omitnan');
Spectral_Flux.KE_Flux_fd_flip = mean(KE_Flux_fd_flip,3, 'omitnan');
Spectral_Flux.qq_Flux_fd_flip = mean(qq_Flux_fd_flip,3, 'omitnan');
Spectral_Flux.vortvort_Flux_fd_flip = mean(vortvort_Flux_fd_flip,3, 'omitnan');


Spectral_Flux.KE_Flux_snapshots = KE_Flux;
Spectral_Flux.qq_Flux_snapshots = qq_Flux;
Spectral_Flux.vortvort_Flux_snapshots = vortvort_Flux;

Spectral_Flux.Wavenumber = kr(:,1,1);

%% Save structure functions ERRORS across snapshots (standard deviation)
Spectral_Flux.KE_Flux_std = std(KE_Flux,0,3, 'omitnan');
Spectral_Flux.qq_Flux_std = std(qq_Flux,0,3, 'omitnan');
Spectral_Flux.vortvort_Flux_std = std(vortvort_Flux,0,3, 'omitnan');


%% Save PE and KE diagnostics to structure

% Spectral_Flux.qq_Work_Rate = qq_Work_Rate;
% Spectral_Flux.qq_Drag_Rate = qq_Drag_Rate;
% Spectral_Flux.qq_Diss_Rate = qq_Diss_Rate;

Spectral_Flux.Time_of_snapshots = time;

%% Save Structure Functions and their errors into a .mat file

save("Spectral_Fluxes_" + experiment_name +"_256.mat",'-struct','Spectral_Flux','-v7.3');

%% Plot spectral fluxes in spectral space

H1_over_H = 0.2/1;
H2_over_H = 0.8/1;

figure
subplot(3,2,1); imagesc(H1_over_H*S.q(:,:,1)); colorbar;
subplot(3,2,2); imagesc(H2_over_H*S.q(:,:,2)); colorbar;
subplot(3,2,3); imagesc(S.psi(:,:,1)); colorbar;
subplot(3,2,4); imagesc(S.psi(:,:,2)); colorbar;
subplot(3,2,5); imagesc(H1_over_H*(S.u(:,:,1).^2 + S.v(:,:,1).^2)); colorbar;
subplot(3,2,6); imagesc(H2_over_H*(S.u(:,:,2).^2 + S.v(:,:,2).^2)); colorbar;

%% Validative analysis of new saved spectral fluxes

figure
subplot(1,2,1)
semilogx(Spectral_Flux.Wavenumber,H1_over_H*Spectral_Flux.qq_Flux(:,1),'k-')
hold on
semilogx(Spectral_Flux.Wavenumber,H2_over_H*Spectral_Flux.qq_Flux(:,2),'r-')
semilogx(Spectral_Flux.Wavenumber, ...
    H1_over_H*Spectral_Flux.qq_Flux(:,1)+H1_over_H*Spectral_Flux.qq_Flux_std(:,1),'k:')
semilogx(Spectral_Flux.Wavenumber, ...
    H1_over_H*Spectral_Flux.qq_Flux(:,1)-H1_over_H*Spectral_Flux.qq_Flux_std(:,1),'k:')
semilogx(Spectral_Flux.Wavenumber, ...
    H2_over_H*Spectral_Flux.qq_Flux(:,2)+H2_over_H*Spectral_Flux.qq_Flux_std(:,2),'r:')
semilogx(Spectral_Flux.Wavenumber, ...
    H2_over_H*Spectral_Flux.qq_Flux(:,2)-H2_over_H*Spectral_Flux.qq_Flux_std(:,2),'r:')
% semilogx([7 max(Spectral_Flux.Wavenumber)], ...
%     Spectral_Flux.qq_Diss_Rate*[1 1], 'r--')
% semilogx([7 max(Spectral_Flux.Wavenumber)], ...
%     mean(Spectral_Flux.qq_Diss_Rate, 'omitnan')*[1 1], 'r-','LineWidth',4)
% semilogx([min(Spectral_Flux.Wavenumber) 7], ...
%     -Spectral_Flux.qq_Drag_Rate*[1 1], 'r--')
% semilogx([7 max(Spectral_Flux.Wavenumber(:,1,1))], ...
%     5e-6*[1 1], 'b-','LineWidth',4)
% semilogx([7 max(Spectral_Flux.Wavenumber)], ...
%     (Spectral_Flux.qq_Diss_Rate + Spectral_Flux.qq_Drag_Rate)*[1 1], 'b--')
% semilogx([7 7], ...
%     [-Spectral_Flux.qq_Drag_Rate Spectral_Flux.qq_Diss_Rate], 'r--')
title('Spectral Flux of Vorticity Variance (2$$\eta$$=$$\zeta^2$$)', Interpreter='latex')
xlim([min(Spectral_Flux.Wavenumber(:,1)) max(Spectral_Flux.Wavenumber(:,1))])


subplot(1,2,2)
top_layer = semilogx(Spectral_Flux.Wavenumber,H1_over_H*Spectral_Flux.KE_Flux(:,1),'k-');
hold on
bottom_layer = semilogx(Spectral_Flux.Wavenumber,H2_over_H*Spectral_Flux.KE_Flux(:,2),'r-');
semilogx(Spectral_Flux.Wavenumber, ...
    H1_over_H*Spectral_Flux.KE_Flux(:,1)+H1_over_H*Spectral_Flux.KE_Flux_std(:,1),'k:')
semilogx(Spectral_Flux.Wavenumber, ...
    H1_over_H*Spectral_Flux.KE_Flux(:,1)-H1_over_H*Spectral_Flux.KE_Flux_std(:,1),'k:')
semilogx(Spectral_Flux.Wavenumber, ...
    H2_over_H*Spectral_Flux.KE_Flux(:,2)+H2_over_H*Spectral_Flux.KE_Flux_std(:,2),'r:')
semilogx(Spectral_Flux.Wavenumber, ...
    H2_over_H*Spectral_Flux.KE_Flux(:,2)-H2_over_H*Spectral_Flux.KE_Flux_std(:,2),'r:')
% semilogx([7 max(Spectral_Flux.Wavenumber)], ...
%     0.5*Spectral_Flux.qq_Diss_Rate*[1 1], 'r--')
% semilogx([min(Spectral_Flux.Wavenumber) 14], ...
%     -0.5*Spectral_Flux.qq_Drag_Rate*[1 1], 'r--')
% semilogx([7 max(Spectral_Flux.Wavenumber)], ...
%     5e-6*[1 1], 'b-','LineWidth',4)
title('Spectral Flux of Kinetic Energy', Interpreter='latex')
xlim([min(Spectral_Flux.Wavenumber(:,1)) max(Spectral_Flux.Wavenumber(:,1))])
legend([top_layer, bottom_layer], 'Top layer', 'Bottom layer')


%% Plot concise spectral fluxes in spectral space

figure
subplot(1,2,1)
semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux,'k-','LineWidth',2)
hold on
semilogx(Spectral_Flux.Wavenumber, ...
    Spectral_Flux.qq_Flux+Spectral_Flux.qq_Flux_std,'k:')
semilogx(Spectral_Flux.Wavenumber, ...
    Spectral_Flux.qq_Flux-Spectral_Flux.qq_Flux_std,'k:')
semilogx(Spectral_Flux.Wavenumber, ...
    quantile(squeeze(Spectral_Flux.qq_Flux_snapshots(:,1,:)),0.75,2),'b--')
semilogx(Spectral_Flux.Wavenumber, ...
    quantile(squeeze(Spectral_Flux.qq_Flux_snapshots(:,1,:)),0.25,2),'b--')
% semilogx([7 max(Spectral_Flux.Wavenumber)], ...
%     mean(Spectral_Flux.qq_Diss_Rate, 'omitnan')*[1 1], 'r-','LineWidth',4)
% semilogx([min(Spectral_Flux.Wavenumber) 7], ...
%     -mean(Spectral_Flux.qq_Drag_Rate, 'omitnan')*[1 1], 'r-','LineWidth',4)
semilogx([min(Spectral_Flux.Wavenumber) max(Spectral_Flux.Wavenumber)], ...
       [0 0], 'k-','LineWidth',1)
% semilogx([7 7], ...
%     [-mean(Spectral_Flux.qq_Drag_Rate, 'omitnan') mean(Spectral_Flux.qq_Diss_Rate, 'omitnan')], 'r--')
title('Flux of buoyancy variance')
xlim([min(Spectral_Flux.Wavenumber) max(Spectral_Flux.Wavenumber)])



subplot(1,2,2)
semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux,'k-','LineWidth',2)
hold on
semilogx(Spectral_Flux.Wavenumber, ...
    Spectral_Flux.KE_Flux+Spectral_Flux.KE_Flux_std,'k:')
semilogx(Spectral_Flux.Wavenumber, ...
    Spectral_Flux.KE_Flux-Spectral_Flux.KE_Flux_std,'k:')
semilogx([min(Spectral_Flux.Wavenumber) max(Spectral_Flux.Wavenumber)], ...
       [0 0], 'k-','LineWidth',1)
% semilogx([min(Spectral_Flux.Wavenumber) 7], ...
%     -0.5*mean(Spectral_Flux.qq_Drag_Rate, 'omitnan')*[1 1], 'r-','LineWidth',4)
% semilogx([7 max(Spectral_Flux.Wavenumber)], ...
%     0.5*mean(Spectral_Flux.qq_Diss_Rate, 'omitnan')*[1 1], 'r-','LineWidth',4)
title('Kinetic Energy Flux')
xlim([min(Spectral_Flux.Wavenumber) max(Spectral_Flux.Wavenumber)])


%% Plot comparison of spectral and finite difference methods for calculating gradients

alpha = 0.9;

figure
subplot(1,2,1)
% First plot the dissipation rates of buoyancy variance at large and small
% scales
% semilogx([7 max(Spectral_Flux.Wavenumber)], ...
%     mean(Spectral_Flux.qq_Diss_Rate, 'omitnan')*[1 1], 'color',[0,0,0]+alpha,'LineWidth',4)
hold on
% semilogx([min(Spectral_Flux.Wavenumber) 7], ...
%     -mean(Spectral_Flux.qq_Drag_Rate, 'omitnan')*[1 1], 'color',[0,0,0]+alpha,'LineWidth',4)
semilogx([min(Spectral_Flux.Wavenumber) max(Spectral_Flux.Wavenumber)], ...
       [0 0], 'k-','LineWidth',1)
% semilogx([7 7], ...
%     [-mean(Spectral_Flux.qq_Drag_Rate, 'omitnan') mean(Spectral_Flux.qq_Diss_Rate, 'omitnan')],  'color',[0,0,0]+alpha)

% Plot the spectral fluxes
a1=semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux,'k-','LineWidth',4);
a2=semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux_fd,'b-','LineWidth',2);
a3=semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.qq_Flux_fd_flip,'r-','LineWidth',2);
title('Flux of buoyancy variance')
xlim([min(Spectral_Flux.Wavenumber) max(Spectral_Flux.Wavenumber)])
%legend([a1 a2 a3], 'Spectral derivatives', 'Finite diff. (small-k zeroed)', 'Finite diff. (large-k zeroed)')
set(gca,'FontSize', 18)

subplot(1,2,2)
% First plot the dissipation rates of buoyancy variance at large and small
% scales
% semilogx([7 max(Spectral_Flux.Wavenumber)], ...
%     0.5*mean(Spectral_Flux.qq_Diss_Rate, 'omitnan')*[1 1], 'color',[0,0,0]+alpha,'LineWidth',4)
hold on
semilogx([min(Spectral_Flux.Wavenumber) 7], ...
    -0.5*mean(Spectral_Flux.qq_Drag_Rate, 'omitnan')*[1 1], 'color',[0,0,0]+alpha,'LineWidth',4)
semilogx([min(Spectral_Flux.Wavenumber) max(Spectral_Flux.Wavenumber)], ...
       [0 0], 'k-','LineWidth',1)
% semilogx([7 7], ...
%     0.5*[-mean(Spectral_Flux.qq_Drag_Rate, 'omitnan') mean(Spectral_Flux.qq_Diss_Rate, 'omitnan')],  'color',[0,0,0]+alpha)

% Plot the spectral fluxes
a1=semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux,'k-','LineWidth',4);
a2=semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux_fd,'b-','LineWidth',2);
a3=semilogx(Spectral_Flux.Wavenumber,Spectral_Flux.KE_Flux_fd_flip,'r-','LineWidth',2);
title('Flux of kinetic energy variance')
xlim([min(Spectral_Flux.Wavenumber) max(Spectral_Flux.Wavenumber)])
set(gca,'FontSize', 18)




