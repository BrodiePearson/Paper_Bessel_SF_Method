%% Calculation of structure functions from snapshots of 2D turbulence
% These structure functions are saved as a structure in a .mat file

%% Load data from 2D simulations
clear all

% Define the beta values for the set of 4 experiments
% Strong meridional beta (SB)
% Strong diagonal beta (SBd)
% Weak meridional beta (WB
% Weak diagonal beta (WBd)
% No beta (NB)
experiment = ["Beta_0.0","Beta_0.5"];
beta = ["0.0", "0.5"];
%betax = ["0.0", "7.071067811865475", "0.0", "0.7071067811865475", "0.0"];

for exp_no = 1:size(experiment,2)
    
    experiment_name = "SurfaceQG_n_512_visc_1.0e-17_order_8_hypovisc_0.005_"+ ...
        "order_-2_kf_7.0_F_1.0e-5_endtime_30000.0_beta_"+beta(exp_no);
    filedir = "../simulations/Output/Data/" + experiment_name + "/";
    filename = filedir + experiment_name + "_*.mat";
    files = dir(filedir + experiment_name + "_*.mat");
    
    % Sort the snapshots from earliest to latest
    temp1 = struct2table(files); % convert the struct array to a table
    temp2 = sortrows(temp1, 'datenum'); % sort the table by 'DOB'
    files = table2struct(temp2);
    
    clear temp1 temp2
    
    file_no = length(files);
    
    direction = ["Zonal", "Meridional", "Diagonal", "Off-Diagonal"];
    
    for dir_no = 1:size(direction,2)
        
        %% Load a data file to set size of structure function arrays before entering
        %  parallel loop
        
        S1 = load(filedir + files(1).name);
        b = S1.b;
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
        
        for nn=1:file_no
            
            tmp_SF_bb=zeros(max_offset,1);
            tmp_SF_LL=zeros(max_offset,1);
            tmp_SF_TT=zeros(max_offset,1);
            tmp_SF_LLL=zeros(max_offset,1);
            tmp_SF_TTL=zeros(max_offset,1);
            tmp_SF_bbL=zeros(max_offset,1);
            tmp_SF_bbadv=zeros(max_offset,1);
            tmp_SF_uuadv=zeros(max_offset,1);
            tmp_SF_vvadv=zeros(max_offset,1);
            tmp_DR_BIN_COUNT=zeros(size(DR_BINS'));
            
            S = load(filedir + files(nn).name)
            b = S.b;
            x = [S.dx/2:S.dx:S.Lx-S.dx/2];
            y = [S.dy/2:S.dy:S.Lx-S.dy/2];
            u = S.u;
            v = S.v;
            buoyancy_work = S.buoyancy_work;
            buoyancy_dissipation_hypoviscosity = S.buoyancy_dissipation_hypoviscosity;
            buoyancy_dissipation_hyperviscosity = S.buoyancy_dissipation_hyperviscosity;
            KE = S.kinetic_energy;
            % Tile x and y positions onto grid
            x=repmat(x',1,size(x,2));
            y=repmat(y,size(y,2),1);
            
            
            
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
            dbdx=(circshift(b,-1,1)-circshift(b,1,1))./(2*x_sep);
            dbdy=(circshift(b,-1,2)-circshift(b,1,2))./(2*y_sep);
            
            % Calculate buoyancy and velocity advection terms
            uadv=u.*dudx+v.*dudy;
            vadv=u.*dvdx+v.*dvdy;
            badv=u.*dbdx+v.*dbdy;
            
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
                u_diff=(circshift(u,[-ii,-jj])-u);
                v_diff=(circshift(v,[-ii,-jj])-v);
                b_diff=(circshift(b,[-ii,-jj])-b);
                uL_diff=(circshift(uL,[-ii,-jj])-uL);
                uT_diff=(circshift(uT,[-ii,-jj])-uT);
                uadv_diff=(circshift(uadv,[-ii,-jj])-uadv);
                vadv_diff=(circshift(vadv,[-ii,-jj])-vadv);
                badv_diff=(circshift(badv,[-ii,-jj])-badv);
                
                tmp_SF_bb(k)=tmp_SF_bb(k)+nanmean(b_diff(:).^2);
                tmp_SF_uuadv(k)=tmp_SF_uuadv(k)+nanmean(u_diff(:).*uadv_diff(:));
                tmp_SF_vvadv(k)=tmp_SF_vvadv(k)+nanmean(v_diff(:).*vadv_diff(:));
                tmp_SF_bbadv(k)=tmp_SF_bbadv(k)+nanmean(b_diff(:).*badv_diff(:));
                tmp_SF_LLL(k)=tmp_SF_LLL(k)+nanmean(uL_diff(:).^3);
                tmp_SF_TTL(k)=tmp_SF_TTL(k)+nanmean(uL_diff(:).*uT_diff(:).^2);
                tmp_SF_bbL(k)=tmp_SF_bbL(k)+nanmean(uL_diff(:).*b_diff(:).^2);
                tmp_SF_LL(k)=tmp_SF_LL(k)+nanmean(uL_diff(:).^2);
                tmp_SF_TT(k)=tmp_SF_TT(k)+nanmean(uT_diff(:).^2);
                
                tmp_DR_BIN_COUNT(k)=tmp_DR_BIN_COUNT(k)+1;
            end
            
            DR_bin_count(:,nn) = tmp_DR_BIN_COUNT;
            SF_bb(:,nn) = tmp_SF_bb./tmp_DR_BIN_COUNT;
            SF_uuadv(:,nn) = tmp_SF_uuadv./tmp_DR_BIN_COUNT;
            SF_vvadv(:,nn) = tmp_SF_vvadv./tmp_DR_BIN_COUNT;
            SF_bbadv(:,nn) = tmp_SF_bbadv./tmp_DR_BIN_COUNT;
            SF_LLL(:,nn) = tmp_SF_LLL./tmp_DR_BIN_COUNT;
            SF_TTL(:,nn) = tmp_SF_TTL./tmp_DR_BIN_COUNT;
            SF_bbL(:,nn) = tmp_SF_bbL./tmp_DR_BIN_COUNT;
            SF_LL(:,nn) = tmp_SF_LL./tmp_DR_BIN_COUNT;
            SF_TT(:,nn) = tmp_SF_TT./tmp_DR_BIN_COUNT;
            
            bb_Work_Rate(nn,1) = buoyancy_work;
            bb_Drag_Rate(nn,1) = buoyancy_dissipation_hypoviscosity;
            bb_Diss_Rate(nn,1) = buoyancy_dissipation_hyperviscosity;
            Kinetic_Energy(nn,1) = KE;
            time(nn)=files(nn).datenum;
        end
        
        %% Save structure functions averaged across snapshots to a structure
        % Also save structure functions for each snapshot for uncert. diagnosis
        
        Structure_Function.bb = nanmean(SF_bb,2);
        Structure_Function.uuadv = nanmean(SF_uuadv,2);
        Structure_Function.vvadv = nanmean(SF_vvadv,2);
        Structure_Function.bbadv = nanmean(SF_bbadv,2);
        Structure_Function.LLL = nanmean(SF_LLL,2);
        Structure_Function.TTL = nanmean(SF_TTL,2);
        Structure_Function.bbL = nanmean(SF_bbL,2);
        Structure_Function.LL = nanmean(SF_LL,2);
        Structure_Function.TT = nanmean(SF_TT,2);
        Structure_Function.R = DR_BINS';
        Structure_Function.R_DIAGONAL = DR_BINS_DIAGONAL';
        
        Structure_Function.bb_snapshots = SF_bb;
        Structure_Function.uuadv_snapshots = SF_uuadv;
        Structure_Function.vvadv_snapshots = SF_vvadv;
        Structure_Function.bbadv_snapshots = SF_bbadv;
        Structure_Function.LLL_snapshots = SF_LLL;
        Structure_Function.TTL_snapshots = SF_TTL;
        Structure_Function.bbL_snapshots = SF_bbL;
        Structure_Function.LL_snapshots = SF_LL;
        Structure_Function.TT_snapshots = SF_TT;
        Structure_Function.Time_of_snapshots = time;
        
        %% Save enstrophy and energy diagnostics to structure
        
        Structure_Function.bb_Work_Rate = bb_Work_Rate;
        Structure_Function.bb_Drag_Rate = bb_Drag_Rate;
        Structure_Function.bb_Diss_Rate = bb_Diss_Rate;
        
        Structure_Function.Kinetic_Energy = Kinetic_Energy;
        
        
        %% Save Structure Functions and their errors into a .mat file
        
        save("Structure_Functions_SQG_"+experiment(exp_no)+"_"+direction(dir_no)+".mat",'-struct','Structure_Function','-v7.3');
        
    end
end
