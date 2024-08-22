%% Calculation of structure functions from snapshots of 2D turbulence
% These structure functions are saved as a structure in a .mat file

%% Load data from 2D simulations
clear all

% Define the beta values for the set of 4 experiments
% Strong meridional beta (SB)
% Strong diagonal beta (SBd)
% Weak meridional beta (WB)
% Weak diagonal beta (WBd)
% No beta (NB)
experiment = ["SB"]; %["SB","SBd","WB","WBd","NB"];
betay = ["10.0"]; %["10.0", "7.071067811865475", "1.0", "0.7071067811865475", "0.0"];
betax = ["0.0"]; %["0.0", "7.071067811865475", "0.0", "0.7071067811865475", "0.0"];

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

direction = ["Zonal", "Meridional", "Diagonal"]; %["Zonal", "Meridional", "Diagonal", "Off-Diagonal"];

for dir_no = 1:size(direction,2)

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

clear S1 x y q u v

%% Calculate structure functions for each snapshot using parallel loops

parfor nn=1:file_no
    
    tmp_SF_vortvort=zeros(max_offset,1);
    tmp_SF_LL=zeros(max_offset,1);
    tmp_SF_TT=zeros(max_offset,1);
    tmp_SF_LLL=zeros(max_offset,1);
    tmp_SF_TTL=zeros(max_offset,1);
    tmp_SF_TTT=zeros(max_offset,1);
    tmp_SF_vortvortL=zeros(max_offset,1);
    tmp_SF_vortvort=zeros(max_offset,1);
    tmp_SF_qqadv=zeros(max_offset,1);
    tmp_SF_uuadv=zeros(max_offset,1);
    tmp_SF_vvadv=zeros(max_offset,1);
    tmp_SF_uuadv_ududx=zeros(max_offset,1);
    tmp_SF_vvadv_udvdx=zeros(max_offset,1);
    tmp_SF_uuadv_vdudy=zeros(max_offset,1);
    tmp_SF_vvadv_vdvdy=zeros(max_offset,1);
    tmp_SF_uuadv_duvdy=zeros(max_offset,1);
    tmp_SF_uuadv_udiv=zeros(max_offset,1);
    tmp_SF_vvadv_duvdx=zeros(max_offset,1);
    tmp_SF_vvadv_vdiv=zeros(max_offset,1);
    tmp_SF_vortvortadv=zeros(max_offset,1);
    tmp_DR_BIN_COUNT=zeros(size(DR_BINS'));
    
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
    duvdx=(circshift(u.*v,-1,1)-circshift(u.*v,1,1))./(2*x_sep);
    duvdy=(circshift(u.*v,-1,2)-circshift(u.*v,1,2))./(2*y_sep);
    
    % Calculate vorticity, vorticity gradiants, and advection terms
    vort=dvdx-dudy;
    divergence = dudx + dvdy;
    dvortdx=(circshift(vort,-1,1)-circshift(vort,1,1))./(2*x_sep);
    dvortdy=(circshift(vort,-1,2)-circshift(vort,1,2))./(2*y_sep);
    uadv_ududx = u.*dudx;
    uadv_vdudy = v.*dudy;
    vadv_udvdx = u.*dvdx;
    vadv_vdvdy = v.*dvdy;
    uadv_udiv = u.*divergence;
    vadv_vdiv = v.*divergence;
    uadv=u.*dudx+v.*dudy;
    vadv=u.*dvdx+v.*dvdy;
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
        u_diff=(circshift(u,[-ii,-jj])-u);
        v_diff=(circshift(v,[-ii,-jj])-v);
        vort_diff=(circshift(vort,[-ii,-jj])-vort);
        uL_diff=(circshift(uL,[-ii,-jj])-uL);
        uT_diff=(circshift(uT,[-ii,-jj])-uT);
        uadv_diff=(circshift(uadv,[-ii,-jj])-uadv);
        vadv_diff=(circshift(vadv,[-ii,-jj])-vadv);
        uadv_ududx_diff=(circshift(uadv_ududx,[-ii,-jj])-uadv_ududx);
        uadv_vdudy_diff=(circshift(uadv_vdudy,[-ii,-jj])-uadv_vdudy);
        uadv_udiv_diff=(circshift(uadv_udiv,[-ii,-jj])-uadv_udiv);
        vadv_udvdx_diff=(circshift(vadv_udvdx,[-ii,-jj])-vadv_udvdx);
        vadv_vdvdy_diff=(circshift(vadv_vdvdy,[-ii,-jj])-vadv_vdvdy);
        vadv_vdiv_diff=(circshift(vadv_vdiv,[-ii,-jj])-vadv_vdiv);
        uadv_duvdy_diff = (circshift(duvdy,[-ii,-jj])-duvdy);
        vadv_duvdx_diff = (circshift(duvdx,[-ii,-jj])-duvdx);
        vortadv_diff=(circshift(vortadv,[-ii,-jj])-vortadv);
        
        tmp_SF_vortvort(k)=tmp_SF_vortvort(k)+mean(vort_diff(:).^2,"omitnan");
        tmp_SF_uuadv(k)=tmp_SF_uuadv(k)+mean(u_diff(:).*uadv_diff(:),"omitnan");
        tmp_SF_vvadv(k)=tmp_SF_vvadv(k)+mean(v_diff(:).*vadv_diff(:),"omitnan");
        tmp_SF_uuadv_ududx(k)=tmp_SF_uuadv_ududx(k)+mean(u_diff(:).*uadv_ududx_diff(:),"omitnan");
        tmp_SF_uuadv_vdudy(k)=tmp_SF_uuadv_vdudy(k)+mean(u_diff(:).*uadv_vdudy_diff(:),"omitnan");
        tmp_SF_vvadv_udvdx(k)=tmp_SF_vvadv_udvdx(k)+mean(v_diff(:).*vadv_udvdx_diff(:),"omitnan");
        tmp_SF_vvadv_vdvdy(k)=tmp_SF_vvadv_vdvdy(k)+mean(v_diff(:).*vadv_vdvdy_diff(:),"omitnan");
        tmp_SF_uuadv_duvdy(k)=tmp_SF_uuadv_duvdy(k)+mean(u_diff(:).*uadv_duvdy_diff(:),"omitnan");
        tmp_SF_vvadv_duvdx(k)=tmp_SF_vvadv_duvdx(k)+mean(v_diff(:).*vadv_duvdx_diff(:),"omitnan");
        tmp_SF_uuadv_udiv(k)=tmp_SF_uuadv_udiv(k)+mean(u_diff(:).*uadv_udiv_diff(:),"omitnan");
        tmp_SF_vvadv_vdiv(k)=tmp_SF_vvadv_vdiv(k)+mean(v_diff(:).*vadv_vdiv_diff(:),"omitnan");
        tmp_SF_vortvortadv(k)=tmp_SF_vortvortadv(k)+mean(vort_diff(:).*vortadv_diff(:),"omitnan");
        tmp_SF_LLL(k)=tmp_SF_LLL(k)+mean(uL_diff(:).^3,"omitnan");
        tmp_SF_TTL(k)=tmp_SF_TTL(k)+mean(uL_diff(:).*uT_diff(:).^2,"omitnan");
        tmp_SF_TTT(k)=tmp_SF_TTT(k)+mean(uT_diff(:).^3,"omitnan");
        tmp_SF_vortvortL(k)=tmp_SF_vortvortL(k)+mean(uL_diff(:).*vort_diff(:).^2,"omitnan");
        tmp_SF_LL(k)=tmp_SF_LL(k)+mean(uL_diff(:).^2,"omitnan");
        tmp_SF_TT(k)=tmp_SF_TT(k)+mean(uT_diff(:).^2,"omitnan");
        
        tmp_DR_BIN_COUNT(k)=tmp_DR_BIN_COUNT(k)+1;
    end
    
    DR_bin_count(:,nn) = tmp_DR_BIN_COUNT;
    SF_vortvort(:,nn) = tmp_SF_vortvort./tmp_DR_BIN_COUNT;
    SF_uuadv(:,nn) = tmp_SF_uuadv./tmp_DR_BIN_COUNT;
    SF_vvadv(:,nn) = tmp_SF_vvadv./tmp_DR_BIN_COUNT;
    SF_uuadv_ududx(:,nn) = tmp_SF_uuadv_ududx./tmp_DR_BIN_COUNT;
    SF_uuadv_vdudy(:,nn) = tmp_SF_uuadv_vdudy./tmp_DR_BIN_COUNT;
    SF_vvadv_udvdx(:,nn) = tmp_SF_vvadv_udvdx./tmp_DR_BIN_COUNT;
    SF_vvadv_vdvdy(:,nn) = tmp_SF_vvadv_vdvdy./tmp_DR_BIN_COUNT;
    SF_uuadv_duvdy(:,nn) = tmp_SF_uuadv_duvdy./tmp_DR_BIN_COUNT;
    SF_vvadv_duvdx(:,nn) = tmp_SF_vvadv_duvdx./tmp_DR_BIN_COUNT;
    SF_uuadv_udiv(:,nn) = tmp_SF_uuadv_udiv./tmp_DR_BIN_COUNT;
    SF_vvadv_vdiv(:,nn) = tmp_SF_vvadv_vdiv./tmp_DR_BIN_COUNT;
    SF_vortvortadv(:,nn) = tmp_SF_vortvortadv./tmp_DR_BIN_COUNT;
    SF_LLL(:,nn) = tmp_SF_LLL./tmp_DR_BIN_COUNT;
    SF_TTL(:,nn) = tmp_SF_TTL./tmp_DR_BIN_COUNT;
    SF_TTT(:,nn) = tmp_SF_TTT./tmp_DR_BIN_COUNT;
    SF_vortvortL(:,nn) = tmp_SF_vortvortL./tmp_DR_BIN_COUNT;
    SF_LL(:,nn) = tmp_SF_LL./tmp_DR_BIN_COUNT;
    SF_TT(:,nn) = tmp_SF_TT./tmp_DR_BIN_COUNT;
    
    Enstrophy_Work_Rate(nn,1) = Work_Rate_ENS; 
    Enstrophy_Drag_Rate(nn,1) = Drag_Rate_ENS; 
    Enstrophy_Diss_Rate(nn,1) = Diss_Rate_ENS; 
    Energy_Work_Rate(nn,1) = Work_Rate_KE; 
    Energy_Drag_Rate(nn,1) = Drag_Rate_KE; 
    Energy_Diss_Rate(nn,1) = Diss_Rate_KE;
    Energy_Diagnostic(nn,1) = Energy;
    Energy_Derived(nn,1) = mean(0.5*(u(:).*u(:) + v(:).*v(:)),"omitnan");
    Enstrophy_Derived(nn,1) = mean(0.5*(q(:).*q(:)),"omitnan");
    
    time(nn)=files(nn).datenum;    
end

%% Save structure functions averaged across snapshots to a structure
% Also save structure functions for each snapshot for uncert. diagnosis

Structure_Function.vortvort = mean(SF_vortvort,2,"omitnan");
Structure_Function.uuadv = mean(SF_uuadv,2,"omitnan");
Structure_Function.vvadv = mean(SF_vvadv,2,"omitnan");
Structure_Function.uuadv_ududx = mean(SF_uuadv_ududx,2,"omitnan");
Structure_Function.uuadv_vdudy = mean(SF_uuadv_vdudy,2,"omitnan");
Structure_Function.vvadv_udvdx = mean(SF_vvadv_udvdx,2,"omitnan");
Structure_Function.vvadv_vdvdy = mean(SF_vvadv_vdvdy,2,"omitnan");
Structure_Function.uuadv_duvdy = mean(SF_uuadv_duvdy,2,"omitnan");
Structure_Function.vvadv_duvdx = mean(SF_vvadv_duvdx,2,"omitnan");
Structure_Function.uuadv_udiv = mean(SF_uuadv_udiv,2,"omitnan");
Structure_Function.vvadv_vdiv = mean(SF_vvadv_vdiv,2,"omitnan");
Structure_Function.vortvortadv = mean(SF_vortvortadv,2,"omitnan");
Structure_Function.LLL = mean(SF_LLL,2,"omitnan");
Structure_Function.TTL = mean(SF_TTL,2,"omitnan");
Structure_Function.TTT = mean(SF_TTT,2,"omitnan");
Structure_Function.vortvortL = mean(SF_vortvortL,2,"omitnan");
Structure_Function.LL = mean(SF_LL,2,"omitnan");
Structure_Function.TT = mean(SF_TT,2,"omitnan");
Structure_Function.R = DR_BINS';
Structure_Function.R_DIAGONAL = DR_BINS_DIAGONAL';

Structure_Function.vortvort_snapshots = SF_vortvort;
Structure_Function.uuadv_snapshots = SF_uuadv;
Structure_Function.vvadv_snapshots = SF_vvadv;
Structure_Function.uuadv_ududx_snapshots = SF_uuadv_ududx;
Structure_Function.uuadv_vdudy_snapshots = SF_uuadv_vdudy;
Structure_Function.vvadv_udvdx_snapshots = SF_vvadv_udvdx;
Structure_Function.vvadv_vdvdy_snapshots = SF_vvadv_vdvdy;
Structure_Function.uuadv_duvdy_snapshots = SF_uuadv_duvdy;
Structure_Function.vvadv_duvdx_snapshots = SF_vvadv_duvdx;
Structure_Function.uuadv_udiv_snapshots = SF_uuadv_udiv;
Structure_Function.vvadv_vdiv_snapshots = SF_vvadv_vdiv;
Structure_Function.vortvortadv_snapshots = SF_vortvortadv;
Structure_Function.LLL_snapshots = SF_LLL;
Structure_Function.TTL_snapshots = SF_TTL;
Structure_Function.TTT_snapshots = SF_TTT;
Structure_Function.vortvortL_snapshots = SF_vortvortL;
Structure_Function.LL_snapshots = SF_LL;
Structure_Function.TT_snapshots = SF_TT;
Structure_Function.Time_of_snapshots = time;

%% Save enstrophy and energy diagnostics to structure

Structure_Function.Enstrophy_Work_Rate = Enstrophy_Work_Rate;
Structure_Function.Enstrophy_Drag_Rate = Enstrophy_Drag_Rate;
Structure_Function.Enstrophy_Diss_Rate = Enstrophy_Diss_Rate;
Structure_Function.Enstrophy = Enstrophy_Derived;

Structure_Function.Energy_Work_Rate = Energy_Work_Rate;
Structure_Function.Energy_Drag_Rate = Energy_Drag_Rate;
Structure_Function.Energy_Diss_Rate = Energy_Diss_Rate;
Structure_Function.Energy = Energy_Derived;
Structure_Function.Energy_Diagnostic = Energy_Diagnostic;


%% Save Structure Functions and their errors into a .mat file

save("Structure_Functions_"+experiment(exp_no)+"_"+direction(dir_no)+"with_trackSFs",'-struct','Structure_Function','-v7.3');

end
end



