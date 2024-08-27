%% Code to diagnose SQG structure functions from GeophysicalFlows.jl
clear all

%% User-specified parameters specifying simulation and sampling

nx = "128";
hyperviscosity = "1.0e-17";
hyperviscosity_order = "8";
hypoviscosity = "0.005";
hypoviscosity_order = "-2";
forcing_wavenumber = "7.0";
forcing_rate = "1.0e-5";
endtime = "30000.0";
beta = "5.0"; % Switch beta off for some simulations where it is zero

% Let the model spin up for two eddy turnover times (one eddy turnover time ~500)
% So remove files that are before time = 1000
start_file = 16;
end_file = 38;

direction = ["Zonal", "Meridional", "Diagonal", "Off-Diagonal"];

%% Find data files

experiment_name = "multilayerqg_2layer_512_beta_5.0_Ld_0.35";

filedir = "../simulations/output/data/";

files = dir(filedir + experiment_name +"_*.mat");

% Sort the snapshots from earliest to latest
temp1 = struct2table(files); % convert the struct array to a table
temp2 = sortrows(temp1, 'datenum'); % sort the table by 'DOB'
files = table2struct(temp2);

clear temp1 temp2

files(end_file:end) = [];
files(1:start_file) = [];
file_no = length(files);

for dir_no = 1:size(direction,2)

    %% Load a data file to set size of structure function arrays before entering
    %  parallel loop

    S1 = load(filedir + files(1).name);
    q = S1.q; % buoyancy
    x = [S1.dx/2:S1.dx:S1.Lx-S1.dx/2];
    y = [S1.dy/2:S1.dy:S1.Lx-S1.dy/2];
    u=S1.u;
    v=S1.v;
    % Tile x and y positions onto grid
    x=repmat(x',1,size(x,2));
    y=repmat(y,size(y,2),1);
    clear S1

    % Find max zonal/meridional separation (in indices)
    % (periodicity means the max distance is half the domain size)
    max_offset=size(u,2)/2;

    %Factor of 2 means that first separation points in diagonal are different
    %bin to pure x/y first separations
    r_interval=abs(x(2,1)-x(1,1));
    DR_BINS=(1:max_offset)*r_interval;
    DR_BINS_DIAGONAL=(1:max_offset)*r_interval*sqrt(2); % Diagonal separations are sqrt(2) larger
    DR_BIN_COUNT=zeros(size(DR_BINS));
    r=zeros(size(DR_BINS));

    clear S1 x y b u v

    %% Calculate structure functions for each snapshot using parallel loops
    for layer=1:2
        parfor nn=1:file_no %parfor nn=1:file_no

            tmp_SF_qq=zeros(max_offset,1);
            tmp_SF_vortvort=zeros(max_offset,1);
            tmp_SF_LL=zeros(max_offset,1);
            tmp_SF_TT=zeros(max_offset,1);
            tmp_SF_qqL=zeros(max_offset,1);
            tmp_SF_qqT=zeros(max_offset,1);
            tmp_SF_vortvortL=zeros(max_offset,1);
            tmp_SF_vortvortT=zeros(max_offset,1);
            tmp_SF_LLL=zeros(max_offset,1);
            tmp_SF_TTL=zeros(max_offset,1);
            tmp_SF_LLT=zeros(max_offset,1);
            tmp_SF_qqadv=zeros(max_offset,1);
            tmp_SF_vortvortadv=zeros(max_offset,1);
            tmp_SF_uuadv=zeros(max_offset,1);
            tmp_SF_vvadv=zeros(max_offset,1);
            tmp_DR_BIN_COUNT=zeros(size(DR_BINS'));

            S = load(filedir + files(nn).name)
            q = S.q(:,:,layer); % buoyancy
            x = [S.dx/2:S.dx:S.Lx-S.dx/2];
            y = [S.dy/2:S.dy:S.Lx-S.dy/2];
            u = S.u(:,:,layer);
            v = S.v(:,:,layer);
            %buoyancy_work = S.buoyancy_work;
            %buoyancy_dissipation_hypoviscosity = S.buoyancy_dissipation_hypoviscosity;
            %buoyancy_dissipation_hyperviscosity = S.buoyancy_dissipation_hyperviscosity;
            KE = S.KE(layer);
            % Tile x and y positions onto grid
            x=repmat(x',1,size(x,2));
            y=repmat(y,size(y,2),1);



            %% Calculate gradients and advective terms for new structure functions

            % Calculate gradients in velocity and buoyancy (for advection)
            % Note that for a variable f(i,j), i is the y-component and j is the
            % x-component
            x_sep=abs(x(2,1)-x(1,1));
            y_sep=abs(y(1,2)-y(1,1));
            dudx=(circshift(u,-1,1)-circshift(u,1,1))./(2*x_sep);
            dudy=(circshift(u,-1,2)-circshift(u,1,2))./(2*y_sep);
            dvdx=(circshift(v,-1,1)-circshift(v,1,1))./(2*x_sep);
            dvdy=(circshift(v,-1,2)-circshift(v,1,2))./(2*y_sep);
            dqdx=(circshift(q,-1,1)-circshift(q,1,1))./(2*x_sep);
            dqdy=(circshift(q,-1,2)-circshift(q,1,2))./(2*y_sep);

            % Calculate relative vorticity from velocity fields
            vort=dvdx-dudy;
            dvortdx=(circshift(vort,-1,1)-circshift(vort,1,1))./(2*x_sep);
            dvortdy=(circshift(vort,-1,2)-circshift(vort,1,2))./(2*y_sep);
       
            % Calculate buoyancy and velocity advection terms
            uadv=u.*dudx+v.*dudy;
            vadv=u.*dvdx+v.*dvdy;
            qadv=u.*dqdx+v.*dqdy;
            vortadv=u.*dvortdx+v.*dvortdy;

            %% Calculate Structure functions from data

            for iii=0:max_offset
                if direction(dir_no) == "Meridional"
                    jj=iii;
                    ii=0;
                elseif direction(dir_no) == "Zonal"
                    jj=0;
                    ii=iii;
                elseif direction(dir_no) == "Diagonal"
                    ii=iii;
                    jj=iii;
                elseif direction(dir_no) == "Off-Diagonal"
                    ii=iii;
                    jj=-iii;
                end
                if (ii^2+jj^2)==0
                    continue
                end
                % Diagnose separation distance and angles (for long/transverse)
                k=sqrt(ii^2+jj^2);
                cosine=ii/k;
                sine=jj/k;
                k=iii; % Diagonal components have to have integer numbered bins

                % Convert velocities to longitudinal and transverse
                % UL=U*DX/DR+V*DY/DR and UT=-U*DY/DR+V*DX/DR

                uL=cosine*u+sine*v;
                uT=-sine*u+cosine*v;

                % Calculate differences in variables
                q_diff=(circshift(q,[-ii,-jj])-q);
                vort_diff=(circshift(vort,[-ii,-jj])-vort);
                u_diff=(circshift(u,[-ii,-jj])-u);
                v_diff=(circshift(v,[-ii,-jj])-v);
                uL_diff=(circshift(uL,[-ii,-jj])-uL);
                uT_diff=(circshift(uT,[-ii,-jj])-uT);
                uadv_diff=(circshift(uadv,[-ii,-jj])-uadv);
                vadv_diff=(circshift(vadv,[-ii,-jj])-vadv);
                qadv_diff=(circshift(qadv,[-ii,-jj])-qadv);
                vortadv_diff=(circshift(vortadv,[-ii,-jj])-vortadv);

                tmp_SF_qq(k)=tmp_SF_qq(k)+mean(q_diff(:).^2);
                tmp_SF_vortvort(k)=tmp_SF_vortvort(k)+mean(vort_diff(:).^2);
                tmp_SF_uuadv(k)=tmp_SF_uuadv(k)+mean(u_diff(:).*uadv_diff(:));
                tmp_SF_vvadv(k)=tmp_SF_vvadv(k)+mean(v_diff(:).*vadv_diff(:));
                tmp_SF_qqadv(k)=tmp_SF_qqadv(k)+mean(q_diff(:).*qadv_diff(:));
                tmp_SF_vortvortadv(k)=tmp_SF_vortvortadv(k)+mean(vort_diff(:).*vortadv_diff(:));
                tmp_SF_LLL(k)=tmp_SF_LLL(k)+mean(uL_diff(:).^3);
                tmp_SF_TTL(k)=tmp_SF_TTL(k)+mean(uL_diff(:).*uT_diff(:).^2);
                tmp_SF_LLT(k)=tmp_SF_LLT(k)+mean(uT_diff(:).*uL_diff(:).^2);
                tmp_SF_qqL(k)=tmp_SF_qqL(k)+mean(uL_diff(:).*q_diff(:).^2);
                tmp_SF_qqT(k)=tmp_SF_qqT(k)+mean(uT_diff(:).*q_diff(:).^2);
                tmp_SF_vortvortL(k)=tmp_SF_vortvortL(k)+mean(uL_diff(:).*vort_diff(:).^2);
                tmp_SF_vortvortT(k)=tmp_SF_vortvortT(k)+mean(uT_diff(:).*vort_diff(:).^2);
                tmp_SF_LL(k)=tmp_SF_LL(k)+mean(uL_diff(:).^2);
                tmp_SF_TT(k)=tmp_SF_TT(k)+mean(uT_diff(:).^2);

                tmp_DR_BIN_COUNT(k)=tmp_DR_BIN_COUNT(k)+1;

            end

            DR_bin_count(:,layer,nn) = tmp_DR_BIN_COUNT;
            SF_qq(:,layer,nn) = tmp_SF_qq./tmp_DR_BIN_COUNT;
            SF_vortvort(:,layer,nn) = tmp_SF_vortvort./tmp_DR_BIN_COUNT;
            SF_uuadv(:,layer,nn) = tmp_SF_uuadv./tmp_DR_BIN_COUNT;
            SF_vvadv(:,layer,nn) = tmp_SF_vvadv./tmp_DR_BIN_COUNT;
            SF_qqadv(:,layer,nn) = tmp_SF_qqadv./tmp_DR_BIN_COUNT;
            SF_vortvortadv(:,layer,nn) = tmp_SF_vortvortadv./tmp_DR_BIN_COUNT;
            SF_LLL(:,layer,nn) = tmp_SF_LLL./tmp_DR_BIN_COUNT;
            SF_TTL(:,layer,nn) = tmp_SF_TTL./tmp_DR_BIN_COUNT;
            SF_LLT(:,layer,nn) = tmp_SF_LLT./tmp_DR_BIN_COUNT;
            SF_qqL(:,layer,nn) = tmp_SF_qqL./tmp_DR_BIN_COUNT;
            SF_qqT(:,layer,nn) = tmp_SF_qqT./tmp_DR_BIN_COUNT;
            SF_vortvortL(:,layer,nn) = tmp_SF_vortvortL./tmp_DR_BIN_COUNT;
            SF_vortvortT(:,layer,nn) = tmp_SF_vortvortT./tmp_DR_BIN_COUNT;
            SF_LL(:,layer,nn) = tmp_SF_LL./tmp_DR_BIN_COUNT;
            SF_TT(:,layer,nn) = tmp_SF_TT./tmp_DR_BIN_COUNT;

            %bb_Work_Rate(nn,1) = buoyancy_work;
            %bb_Drag_Rate(nn,1) = buoyancy_dissipation_hypoviscosity;
            %bb_Diss_Rate(nn,1) = buoyancy_dissipation_hyperviscosity;
            Kinetic_Energy(layer,nn) = KE;
            time(nn)=files(nn).datenum;


        end
    end
    %% Save structure functions averaged across snapshots to a structure
    % Also save structure functions for each snapshot for uncert. diagnosis

    Structure_Function.qq = mean(SF_qq,3,'omitnan');
    Structure_Function.vortvort = mean(SF_vortvort,3,'omitnan');
    Structure_Function.uuadv = mean(SF_uuadv,3,'omitnan');
    Structure_Function.vvadv = mean(SF_vvadv,3,'omitnan');
    Structure_Function.qqadv = mean(SF_qqadv,3,'omitnan');
    Structure_Function.vortvortadv = mean(SF_vortvortadv,3,'omitnan');
    Structure_Function.LLL = mean(SF_LLL,3,'omitnan');
    Structure_Function.TTL = mean(SF_TTL,3,'omitnan');
    Structure_Function.LLT = mean(SF_LLT,3,'omitnan');
    Structure_Function.qqL = mean(SF_qqL,3,'omitnan');
    Structure_Function.qqT = mean(SF_qqT,3,'omitnan');
    Structure_Function.vortvortL = mean(SF_vortvortL,3,'omitnan');
    Structure_Function.vortvortT = mean(SF_vortvortT,3,'omitnan');
    Structure_Function.LL = mean(SF_LL,3,'omitnan');
    Structure_Function.TT = mean(SF_TT,3,'omitnan');
    Structure_Function.R = DR_BINS';
    Structure_Function.R_DIAGONAL = DR_BINS_DIAGONAL';

    Structure_Function.qq_snapshots = SF_qq;
    Structure_Function.vortvort_snapshots = SF_vortvort;
    Structure_Function.uuadv_snapshots = SF_uuadv;
    Structure_Function.vvadv_snapshots = SF_vvadv;
    Structure_Function.qqadv_snapshots = SF_qqadv;
    Structure_Function.vortvortadv_snapshots = SF_vortvortadv;
    Structure_Function.LLL_snapshots = SF_LLL;
    Structure_Function.TTL_snapshots = SF_TTL;
    Structure_Function.LLT_snapshots = SF_LLT;
    Structure_Function.qqL_snapshots = SF_qqL;
    Structure_Function.qqT_snapshots = SF_qqT;
    Structure_Function.vortvortL_snapshots = SF_vortvortL;
    Structure_Function.vortvortT_snapshots = SF_vortvortT;
    Structure_Function.LL_snapshots = SF_LL;
    Structure_Function.TT_snapshots = SF_TT;
    Structure_Function.Time_of_snapshots = time;

    %% Save enstrophy and energy diagnostics to structure

    %Structure_Function.bb_Work_Rate = bb_Work_Rate;
    %Structure_Function.bb_Drag_Rate = bb_Drag_Rate;
    %Structure_Function.bb_Diss_Rate = bb_Diss_Rate;

    Structure_Function.Kinetic_Energy = Kinetic_Energy;


    %% Save Structure Functions and their errors into a .mat file

    save("processed_data/Structure_Functions_" + experiment_name +"_"+ direction(dir_no) +".mat",'-struct','Structure_Function','-v7.3');

end
%% Create preliminary plots of structure function data
% First plot the buoyancy and Advection term fields

S = load(filedir + files(file_no).name);
q = S.q; % buoyancy
x = S.dx/2:S.dx:S.Lx-S.dx/2;
y = S.dy/2:S.dy:S.Lx-S.dy/2;
u = S.u;
v = S.v;

figure 

subplot(1,3,1)
imagesc(x(1:end,1),y(1,1:end),q(:,:,1))
colormap gray
colorbar
xlabel('x (m)')
ylabel('y (m)')
title('Vorticity')
set(gca,'fontsize', 17);
ax = gca;
ax.YDir = 'normal';

subplot(1,3,2)
imagesc(x(1:end,1),y(1,1:end),q(:,:,1).^2)
colormap jet
colorbar
xlabel('x (m)')
ylabel('y (m)')
title('Vorticity Variance')
set(gca,'fontsize', 17);
ax = gca;
ax.YDir = 'normal';

subplot(1,3,3)
imagesc(x(1:end,1),y(1,1:end),u(:,:,1).^2+v(:,:,1).^2)
colormap jet
colorbar
xlabel('x (m)')
ylabel('y (m)')
title('Kinetic Energy x 2')
set(gca,'fontsize', 17);
ax = gca;
ax.YDir = 'normal';

%% Calculate divergence of vorticity-vorticity-L structure function

r_interval=DR_BINS(1,2)-DR_BINS(1,1);

dqqLdr=(circshift(Structure_Function.qqL,-1,1)- ...
    circshift(Structure_Function.qqL,1,1))./(2*r_interval);

%% Plot the relation between the new buoyancy SF and the Lindborg 3rd-order SF

R_BINS=DR_BINS';

figure 

subplot(1,3,1)
a1=loglog(R_BINS,-Structure_Function.qqadv,'ko', 'MarkerFaceColor', 'k');
hold on
loglog(R_BINS,Structure_Function.qqadv,'ko')
% loglog([min(Structure_Function.R) 1/7], ...
%     nanmean(Structure_Function.bb_Diss_Rate)*[1 1], 'r-','LineWidth',4)
leg1 = legend(a1, '$\delta q \delta A_{q}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17, 'Location','South');
xlabel('Separation, r (m)')
ylabel('Structure Functions (s^{-3})')
ylim([-1e-5 1e-5])
xlim([1e-2 1e0])
set(gca,'fontsize', 17);
subplot(1,3,3)
a2=loglog(R_BINS,-Structure_Function.qqL./R_BINS,'ko', 'MarkerFaceColor', 'k');
hold on
loglog(R_BINS,Structure_Function.qqL./R_BINS,'ko')
loglog(R_BINS,-Structure_Function.qqT./R_BINS,'ro', 'MarkerFaceColor', 'r')
loglog(R_BINS,Structure_Function.qqT./R_BINS,'ro')
%loglog([min(Structure_Function.R) 1/7], ...
%    nanmean(Structure_Function.qq_Diss_Rate)*[1 1], 'r-','LineWidth',4)
leg1 = legend(a2, '$-(\delta u_L \delta q \delta q)/r$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17, 'Location','South');
xlabel('Separation, r (m)')
ylim([1e-7 1e-5])
xlim([1e-2 1e0])
set(gca,'fontsize', 17);
subplot(1,3,2)
a1=loglog(R_BINS,-0.5*dqqLdr,'ko', 'MarkerFaceColor', 'k');
hold on
loglog(R_BINS,0.5*dqqLdr,'ko')
%loglog([min(Structure_Function.R) 1/7], ...
%    nanmean(Structure_Function.qq_Diss_Rate)*[1 1], 'r-','LineWidth',4)
leg1 = legend(a1, '$-(1/2)\partial(\delta u_L \delta q \delta q)/\partial r$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17, 'Location','South');
xlabel('Separation, r (m)')
ylim([1e-7 1e-5])
xlim([1e-2 1e0])
set(gca,'fontsize', 17);

%% Plot velocity structure functions

figure 

subplot(1,3,1)
a1=loglog(R_BINS,-Structure_Function.uuadv+Structure_Function.vvadv,'ko', 'MarkerFaceColor', 'k');
hold on
loglog(R_BINS,Structure_Function.uuadv+Structure_Function.vvadv,'ko')
%loglog([min(Structure_Function.R) 1/7], ...
%    nanmean(Structure_Function.bb_Diss_Rate)*[1 1], 'r-','LineWidth',4)
leg1 = legend(a1, '$\delta u \cdot \delta A_{u}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17, 'Location','South');
xlabel('Separation, r (m)')
ylabel('Structure Functions (s^{-3})')
ylim([-1e-5 1e-5])
xlim([1e-2 1e0])
set(gca,'fontsize', 17);
subplot(1,3,3)
a2=loglog(R_BINS,-Structure_Function.qqL./R_BINS,'ko', 'MarkerFaceColor', 'k');
hold on
loglog(R_BINS,Structure_Function.qqL./R_BINS,'ko')
loglog(R_BINS,-Structure_Function.qqT./R_BINS,'ro', 'MarkerFaceColor', 'r')
loglog(R_BINS,Structure_Function.qqT./R_BINS,'ro')
%loglog([min(Structure_Function.R) 1/7], ...
%    nanmean(Structure_Function.bb_Diss_Rate)*[1 1], 'r-','LineWidth',4)
leg1 = legend(a2, '$-(\delta u_L \delta q \delta q)/r$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17, 'Location','South');
xlabel('Separation, r (m)')
ylim([1e-7 1e-5])
xlim([1e-2 1e0])
set(gca,'fontsize', 17);
subplot(1,3,2)
a1=loglog(R_BINS,-0.5*dqqLdr,'ko', 'MarkerFaceColor', 'k');
hold on
loglog(R_BINS,0.5*dqqLdr,'ko')
%loglog([min(Structure_Function.R) 1/7], ...
%    nanmean(Structure_Function.bb_Diss_Rate)*[1 1], 'r-','LineWidth',4)
leg1 = legend(a1, '$-(1/2)\partial(\delta u_L \delta q \delta q)/\partial r$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17, 'Location','South');
xlabel('Separation, r (m)')
ylim([1e-7 1e-5])
xlim([1e-2 1e0])
set(gca,'fontsize', 17);


