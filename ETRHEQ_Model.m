% Authors of MATLAB code: Angela Rigden & Guido Salvucci
% Date: August, 2016
%
% This Matlab code is an example of the ETRHEQ method. See the following
% citation for more details on the ETRHEQ method and pre-processing of
% AmeriFlux input/validation data:
%
% Rigden, A. J., and G. D. Salvucci (2016), Stomatal response to decreasing
% humidity implicated in recent decline in U.S. evaporation, Global Change
% Biology, doi:10.1111/gcb.13439.
%
% Please send questions/concerns/suggestions to:
% Angela Rigden at ajrigden@bu.edu
%
% In the following example, the ETRHEQ method is run for one year at one
% AmeriFlux site. The input forcing data is in "sample_data.mat".

% Site name:  Audubon_Research_Ranch, AZ, USA.
% Year:       2006
% Latitude:   31.5907
% Longitude:  -110.5092
% Elevation:  1469 m
% PI:         Tilden Meyers - Tilden.Meyers@noaa.gov - NOAA/ARL
%
% The AmeriFlux data was obtained from http://ameriflux.ornl.gov. Funding
% for AmeriFlux data resources was provided by the U.S. Department of
% Energy's Office of Science.


% ---------- OPTIONS ----------
% Do you want to use observed (1) or modeled (0) downwelling longwave raditation (rld)?
obs_rld = 1; % == 1, use observed rld
             % == 0, use modeled rld
% Do you want to use observed (1) or modeled (0) ground heat flux (ghf)?
obs_ghf = 1; % == 1, use observed ghf
             % == 0, use modeled ghf

% Note that model runs are saved as:
% ['test_run_',num2str(obs_ghf),num2str(obs_rld),'.mat'], i.e.
% 'test_run_00.mat' --> rld and ghf are modeled
% 'test_run_10.mat' --> rld is modeled, ghf is observed
%  etc.
% ----------------------------- 

% Site specific information: Audubon Research Ranch
z_veg    = 0.15;   % z_veg, vegetation height (m)
z_m      = 4;      % z_m, instrument height for humidity and temperature (m)
z_m_wind = 4;      % z_m_wind, instrument height for wind (m)
z_elev   = 1469;   % z_elev, elevation (m)

global r_d cp v k g lv emis eps sb RIMAX

% Assumed constants
kB    = 2;     % kB, kB^-1 ()
emis  = 0.98;  % emis, emissivity of the ground-surface ()
RIMAX = 0.15;  % RIMAX, critical value of flux Richardson number ()
Isoil = 1300;  % Isoil, thermal inertia of soil (J/(m^2 s^(1/2) K))
               % Isoil is only used when modeling ghf (obs_ghf == 0)

% Physical constants
BOLTZMAN = 1.380658e-23;            % BOLTZMAN,  Boltzman constant (J/K)
AVOGADRO = .602214199e24;           % AVOGADRO, Avogadro constant (1/mol)
MD       = 28.9644e-3;              % MD, molar mass dry air (kg/mol)
MV  = 18.0153e-3;                   % MD, molae mass water vapor (kg/mol)
r_v = (AVOGADRO)*(BOLTZMAN) / (MV); % r_v, gas constant for water vapor (J/(kg-K))
r_d = (AVOGADRO)*(BOLTZMAN) / (MD); % r_d, gas constant for dry air (J/(kg-K))
cp  = 7./2*(r_d);                   % cp, specific heat of air (J/(kg-K))
v   = 1.4531e-05;                   % v, kinematic viscosity (m^2/s]
k   = 0.41;                         % k, von Karman constant ()
g   = 9.81;                         % g, gravitational acceleration (m/s^2)
eps = r_d/r_v;                      % eps, ratio of universal gas constant dry air to water vapor ()
lv  = 2.5008e6;                     % lv, latent heat of vaporization (J/kg)
sb  = 5.6704*10^(-8);               % sb, stephan-boltzman constant (W/(m^2 K^4))

% Estimate roughness and displacement heights
z_o  = 0.1*z_veg;     % momentum roughness (m)
d    = 0.7*z_veg;     % displacement hieight (m)
z_ov = (z_o./exp(kB));% roughness for water vapor transport (m)
z_oh = z_ov;          % roughness for heat transfer (m)

% Load half hourly input data, variable: sample_data
load('sample_data.mat');
% sample_data.mat includes three variables:
% sample_data_site_name, indicates AmeriFlux site name
% sample_data_year, indicates year of AmeriFlux data
% sample_data, data (both input and validation)
%       size(sample_data): 48 half hours x 365 days x 14 observed variables
%       The 14 observed variables:
%       1-  t, temperature (K)
%       2-  q, specific humidity (kg/kg)
%       3-  u, windspeed (m/s)
%       4-  p, pressure (Pa)
%       5-  u_str, friction velocity (m/s)
%       6-  le, latent heat flux (W/m^2)
%       7-  sh, sensible heat flux (W/m^2)
%       8-  ghf, ground heat flux (W/m^2)
%       9-  rld, downwelling longwave radiation (W/m^2)
%       10- rlu, upwelling longwave radiation (W/m^2)
%       11- rsd, downwelling solar radiation (W/m^2)
%       12- rsu, upwelling solar radiation (W/m^2)
%       13- rnet, net radiation (W/m^2)
%       14- precip, precipitation (mm)

% Define range of surface resistance to water vapor transport
r_surf = exp(15:-0.25:0); % daily-constant effective surface resistence (s/m)

% Initialize matrices. The size of matrices correspond to
%   48 half hours x 61 values in r_surf x 365 days in year.
%   X_vary_rs store solutions for each value of r_surf.
vertvar_vary_rs = zeros(365,61)+NaN;     % vertvar-, vertical variance of RH averaged over the day ()
ts_vary_rs      = zeros(48,61,365)+NaN;  % ts-, surface temperature (K)
le_vary_rs      = zeros(48,61,365)+NaN;  % le-, latent heat flux (W/m^2)
sh_vary_rs      = zeros(48,61,365)+NaN;  % sh-, sensible heat flux (W/m^2)
rlu_vary_rs     = zeros(48,61,365)+NaN;  % rlu-, upwelling longwave radiation (W/m^2)
u_str_vary_rs   = zeros(48,61,365)+NaN;  % u_str-, friction velocity (m/s)
r_atm_vary_rs   = zeros(48,61,365)+NaN;  % r_atm-, atmospheric resistance (s/m)
L_vary_rs       = zeros(48,61,365)+NaN;  % L-, Obukhov length (m)

% If the day does not have the required input data to run the ETRHEQ
% method, we will put an NaN in the "missing_days" vector, which is
% otherwise filled with ones.
missing_days = zeros(365,1)+1;

% if we are modeling rld/ghf, then the saved modeled results will be useful
% for comparison with observations
rld_etrheq = zeros(48,365)+NaN;
ghf_etrheq = zeros(48,365)+NaN;

% groung heat flux is modeled as a fucntion of air temperature
if obs_ghf == 0
    t = sample_data(:,:,1);
    ghf_mod =  calc_ghf(t,Isoil);
end

for DY = 1:365; % run for all days in year, DY counts days
    
    % Extract data from sample_data for day, DY
    t = sample_data(:,DY,1); % e.g. size(t) = 48 half hours x 1
    q = sample_data(:,DY,2);
    u = sample_data(:,DY,3);
    u(u<=0) = NaN; % u > 0
    p = sample_data(:,DY,4);
    
    rsnet = sample_data(:,DY,11)-...
        sample_data(:,DY,12); % net solar radiation, rsnet (W/m^2)
    rsnet(rsnet<0) = 0; % rsnet >= 0
    
    % Determine if you are using modeled or observed rld/ghf
    if obs_rld == 1      % use observed rld
        rld_tot = sample_data(:,DY,9);
    elseif obs_rld == 0  % use modeled rld
        rld_tot = calc_rld(p,t,q,z_elev);
    end
    % subract out relected longwave radiation
    rld = rld_tot - ((1-emis)*rld_tot);
    
    if obs_ghf == 1      % use observed ghf
        ghf = sample_data(:,DY,8);
    elseif obs_ghf == 0  % use modeled ghf
        ghf = ghf_mod(:,DY);
    end
    
    % put all required ETRHEQ data into one matrix, data
    data = [t,q,u,p,rsnet,rld,ghf];
    
    % if any data is missing for that day, do not run ETRHEQ method for
    % that day (skip)
    if any(isnan(data(:)))
        missing_days(DY) = NaN; % to indicate missing day, replace 1 with NaN
        disp(['... skipping day ',num2str(DY),'/365'])
    else
        
        den = p./(r_d.*t); % density of air, den (kg/m^3)
        t_ref = mean(t);   % reference temperature, t_ref (K)
        
        % If z_m_wind ~= z_m (i.e, the wind is measured at a different
        % height than the q and t, which is often the case at common
        % weather stations), we need to estimate u at z_m. To do this, we
        % assume that the wind at z_m is a fraction of the wind measured at
        % z_m_wind (assuming z_m < z_m_wind).
        
        if z_m == z_m_wind
            num_frac = 8; % number of fractions between 0.1 and 1, num_frac
        else
            num_frac = 8;
        end
        frac_list = linspace(0.1,1,num_frac);
        % if z_m_wind == z_m, frac_list == 1
        % if z_m_wind ~= z_m, frac_list ranges bewteen 0.1 and 1
        
        % Reshape/repmat frac_list to be 48 x 61*num_frac, frac_list_final
        frac_list_final = reshape(permute(repmat(frac_list,[48,1,61]),[1,3,2]),48,61*num_frac);
        
        % Calcuate possible u at z_m values as fractions of u at z_m_wind
        u_pred_zm = bsxfun(@times,frac_list_final,u);
        
        % Repmat r_surf to be 48 x 61*num_frac, r_surf_frac
        r_surf_frac = repmat(ones(48,1)*r_surf,1,num_frac);
        
        % X_frac stores solutions for a range of predicted windspeeds (if
        % z_m ~= z_m_wind), in addition to r_surf values
        % chih- is: log((z_m-d)./z_ov) + stability correction
        [ts_frac, le_frac, sh_frac, rlu_frac, L_frac, chih_frac, ...
            u_str_frac, rib_frac] = calc_ts_bisection(data,...
            r_surf_frac, den, z_m, z_o, z_ov, z_oh, d, u_pred_zm);
        
        % Detemine windspeed at z_m
        % Recall equations:
        % u_zm_wind = (log((z_m_wind-d)./z_o)-PSIMW)*u_str./k
        % u_zm      = (log((z_m-d)./z_o)-PSIM)*u_str./k
        % therefore,
        % u_zm_wind = [(log((z_m_wind-d)./z_o)-PSIMW)./(log((z_m-d)./z_o)-PSIM)]*u_zm
        PSIM  = calc_stability_corr(z_m, d, z_o, L_frac);
        PSIMW = calc_stability_corr(z_m_wind, d, z_o, L_frac);
        frac_new = (log(((z_m_wind-d)./z_o))-PSIMW)./(log((z_m-d)./z_o)-PSIM);
        upred_zm_wind_frac = frac_new.*u_pred_zm;
        upred_zm_wind = reshape(upred_zm_wind_frac,[48,61,num_frac]);
        
        % reshape X_frac to 48 x 61 x num_frac
        ts_frac  = reshape(ts_frac,[48,61,num_frac]);
        le_frac  = reshape(le_frac,[48,61,num_frac]);
        sh_frac  = reshape(sh_frac,[48,61,num_frac]);
        rlu_frac = reshape(rlu_frac,[48,61,num_frac]);
        L_frac   = reshape(L_frac,[48,61,num_frac]);
        chih_frac  = reshape(chih_frac,[48,61,num_frac]);
        u_str_frac = reshape(u_str_frac,[48,61,num_frac]);
        rib_frac   = reshape(rib_frac,[48,61,num_frac]);
        
        ts_final_day = zeros(48,61)+NaN;
        le_final_day = zeros(48,61)+NaN;
        sh_final_day = zeros(48,61)+NaN;
        rlu_final_day = zeros(48,61)+NaN;
        L_final_day = zeros(48,61)+NaN;
        chih_final_day = zeros(48,61)+NaN;
        u_str_final_day = zeros(48,61)+NaN;
        
        % If the critical Ri > RIMAX, set predicted u to NaN;
        upred_zm_wind(rib_frac==RIMAX) = NaN;
        
        for RS = 1:61
            for HR = 1:48
                
                u_hr = u(HR);
                u_pred = squeeze(upred_zm_wind(HR,RS,:));
                
                % Find the frac that best predicts the windspeed at 10 m
                if all(isnan(u_pred))
                    frac_ind = num_frac; % If critical Ri > RIMAX, set to lowest u
                else
                    [~, frac_ind] = min(abs(u_hr-u_pred));
                end
                
                % Save output
                % [48 hours, 61 resistances, num_frac fractions];
                ts_final_day(HR,RS) = ts_frac(HR,RS,frac_ind);
                le_final_day(HR,RS) = le_frac(HR,RS,frac_ind);
                sh_final_day(HR,RS) = sh_frac(HR,RS,frac_ind);
                rlu_final_day(HR,RS)  = rlu_frac(HR,RS,frac_ind);
                L_final_day(HR,RS)    = L_frac(HR,RS,frac_ind);
                chih_final_day(HR,RS) = chih_frac(HR,RS,frac_ind);
                u_str_final_day(HR,RS)= u_str_frac(HR,RS,frac_ind);
            end
        end
        
        
        % Account for condensation by setting r_surf to 1, i.e. r_surf(61), if le<0
        le_thresh = 0; % threshold for condenstaion (W/m^2)
        
        le_con    = le_final_day(:,61)*ones(1,61);
        sh_con    = sh_final_day(:,61)*ones(1,61);
        ts_con    = ts_final_day(:,61)*ones(1,61);
        rlu_con   = rlu_final_day(:,61)*ones(1,61);
        L_con     = L_final_day(:,61)*ones(1,61);
        chih_con  = chih_final_day(:,61)*ones(1,61);
        u_str_con = u_str_final_day(:,61)*ones(1,61);
        
        con = find(le_final_day(:,61) < le_thresh);
        
        le_final_day(con,:)    = le_con(con,:);
        sh_final_day(con,:)    = sh_con(con,:);
        ts_final_day(con,:)    = ts_con(con,:);
        rlu_final_day(con,:)   = rlu_con(con,:);
        L_final_day(con,:)     = L_con(con,:);
        chih_final_day(con,:)  = chih_con(con,:);
        u_str_final_day(con,:) = u_str_con(con,:);
        
        % Set r_surf to inf at night
        sw_thresh = 0.1; % threshold for night (W/m^2)
        
        le_night    = le_final_day(:,1)*ones(1,61);
        h_night     = sh_final_day(:,1)*ones(1,61);
        ts_night    = ts_final_day(:,1)*ones(1,61);
        r_lu_night  = rlu_final_day(:,1)*ones(1,61);
        L_night     = L_final_day(:,1)*ones(1,61);
        chih_night  = chih_final_day(:,1)*ones(1,61);
        u_str_night = u_str_final_day(:,1)*ones(1,61);
        
        r_surf_night = ones(48,1)*r_surf;
        rs_night    = r_surf_night(:,1)*ones(1,61);
        
        swindicator = rsnet*ones(1,61);
        night = find(swindicator < sw_thresh & le_final_day > le_thresh);
        
        le_final_day(night)     = le_night(night);
        sh_final_day(night)     = h_night(night);
        ts_final_day(night)     = ts_night(night);
        rlu_final_day(night)    = r_lu_night(night);
        L_final_day(night)      = L_night(night);
        chih_final_day(night)   = chih_night(night);
        u_str_final_day(night)  = u_str_night(night);
        
        r_surf_night(night) = rs_night(night);
        
        % calculate atmospheric resistance
        r_atm = chih_final_day./(k.*u_str_final_day);
        
        % specify vertical levels to calculate vertical variance over.  In
        % this case, we are using two levels, z1 & z2, which correspond to
        % the surface and the measurement height
        z1 = (z_oh + d);
        z2 = z_m;
        z_levels = cat(3,ones(48,61)*z1, ones(48,61)*z2);
        
        % Estimate RH at z_levels
        RHz = calc_RH_profile(data, u_str_final_day, ts_final_day,...
            sh_final_day, le_final_day, L_final_day, z_levels,...
            r_surf_night, t_ref, den*ones(1,61), ...
            ones(48,61)*z_oh, ones(48,61)*z_ov, d);
        
        % Reshape RH profiles
        RHz = permute(RHz, [3,2,1]);
        
        % Esitmate half_hourly RH profiles at each resistance
        var_RHz = squeeze(nanvar(RHz,0,1))';
        
        % Average the vertical variance over the day (see note below)
        mean_var = nanmean(var_RHz);
        
        % Note: when using two levels, the mean of the 48 variances of the
        % two levels is the same as the sum of squares of the differences,
        % such that the procedure is equivalent to minimizing the
        % root-mean-square of RH(surface)-RH(air). Previous versions (cited
        % below) took variance over twenty equally chi-spaced levels in the
        % boundary layer, but testing revealed that this simpler
        % (two-level) method works just as well.
        
        % Salvucci, G. D., and P. Gentine (2013), Emergent relation between
        % surface vapor conductance and relative humidity profiles yields
        % evaporation rates from weather data, Proceedings of the National
        % Academy of Sciences of the United States of America, 110(16),
        % 6287?6291, doi:10.1073/pnas.1215844110.
        
        % Rigden, A. J., and G. D. Salvucci (2015), Evapotranspiration
        % based on equilibrated relative humidity (ETRHEQ): Evaluation over
        % the continental U.S, Water Resour. Res., 51(4), 2951?2973,
        % doi:10.1002/2014WR016072.
        
        % Store daily results in matrices
        vertvar_vary_rs(DY,:) = mean_var;
        ts_vary_rs(:,:,DY)  = ts_final_day;
        le_vary_rs(:,:,DY)  = le_final_day;
        sh_vary_rs(:,:,DY)  = sh_final_day;
        rlu_vary_rs(:,:,DY) = rlu_final_day;
        u_str_vary_rs(:,:,DY) = u_str_final_day;
        r_atm_vary_rs(:,:,DY) = r_atm;
        L_vary_rs(:,:,DY)     = L_final_day;
        
        % beacuse rld/ghf does not depend on r_surf, but save
        rld_etrheq(:,DY)   = rld;
        
        disp(['... running day ',num2str(DY),'/365'])
        
    end
end

% Calculate precipitation occurance
daily_precip = nansum(sample_data(:,:,14)); % daily precipitaiton, daily_precip (mm)
precip_occ = zeros(size(daily_precip)); % precipitation occurace, precip_occ (0=no precip, 1=precip)
precip_occ(daily_precip>0) = 1;

% ---- main objective function ----
% chooses r_surf (RS indicates index, 1-61) by minimizing the vertical
% variance of RH averaged over the day
RS = windave_var_rs(vertvar_vary_rs', precip_occ);

% X_etrheq are the final ETRHEQ estimate of variable
ts_etrheq = zeros(48,365)+NaN;
le_etrheq = zeros(48,365)+NaN;
sh_etrheq = zeros(48,365)+NaN;
rlu_etrheq   = zeros(48,365)+NaN;
u_str_etrheq = zeros(48,365)+NaN;
r_atm_etrheq = zeros(48,365)+NaN;
L_etrheq     = zeros(48,365)+NaN;

for DY = 1:365
    
    % Diurnal days (48,365)
    if ~isnan(RS(DY))
        ts_etrheq(:,DY)    = ts_vary_rs(:,RS(DY),DY); % all hours, at specific resistance, on specific day, for chosen KB
        le_etrheq(:,DY)    = le_vary_rs(:,RS(DY),DY);
        sh_etrheq(:,DY)    = sh_vary_rs(:,RS(DY),DY);
        rlu_etrheq(:,DY)   = rlu_vary_rs(:,RS(DY),DY);
        u_str_etrheq(:,DY) = u_str_vary_rs(:,RS(DY),DY);
        r_atm_etrheq(:,DY) = r_atm_vary_rs(:,RS(DY),DY);
        L_etrheq(:,DY)     = L_vary_rs(:,RS(DY),DY);
    end
    
end

% X_obs are the observed data at AmeriFlux sites. Recall, these
% are not used in the ETRHEQ method and are only used as validation data.
le_obs = sample_data(:,:,6);
sh_obs = sample_data(:,:,7);
ghf_obs = sample_data(:,:,8);
rld_obs = sample_data(:,:,9);
rlu_obs = sample_data(:,:,10);
rsd_obs = sample_data(:,:,11);
rsu_obs = sample_data(:,:,12);

% X_obs_enbal are the fluxes estimated as a residual of the energy balance,
% e.g. le = rsd - rsu + rld - rlu - sh - ghf
le_obs_enbal = rsd_obs - rsu_obs + rld_obs - rlu_obs - sh_obs - ghf_obs;
sh_obs_enbal = rsd_obs - rsu_obs + rld_obs - rlu_obs - le_obs - ghf_obs;

daily_rmse_le = sqrt(nanmean((nanmean(le_etrheq)-nanmean(le_obs)).^2));
daily_rmse_sh = sqrt(nanmean((nanmean(sh_etrheq)-nanmean(sh_obs)).^2));

save(['test_run_',num2str(obs_ghf),num2str(obs_rld),'.mat'])
% 'test_run_00.mat --> rld and ghf are modeled
% 'test_run_10.mat --> rld is modeled, ghf is observed
%  etc.

% Plot latent heat flux (le)
figure('color','white')
plot(nanmean(le_etrheq),'b-','linewidth',1.5)
hold on
plot(nanmean(le_obs),'r-','linewidth',1.5)
plot(nanmean(le_obs_enbal),'g-','linewidth',1)
ylabel('LE (W/m^2)')
xlabel('Day')
xlim([1 365])
legend('ETRHEQ method','Observed','Observed, energy balance residual','location','northwest')
legend boxoff
title([regexprep(sample_data_site_name,'_',' '), ' (',num2str(sample_data_year),')']);
set(gca,'FontSize',14)

% Plot sensisble heat flux (sh)
figure('color','white')
plot(nanmean(sh_etrheq),'b-','linewidth',1.5)
hold on
plot(nanmean(sh_obs),'r-','linewidth',1.5)
plot(nanmean(sh_obs_enbal),'g-','linewidth',1)
ylabel('SH (W/m^2)')
xlabel('Day')
xlim([1 365])
legend('ETRHEQ method','Observed','Observed, energy balance residual')
legend boxoff
title([regexprep(sample_data_site_name,'_',' '), ' (',num2str(sample_data_year),')']);
set(gca,'FontSize',14)
