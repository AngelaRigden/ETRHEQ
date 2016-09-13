function  [RS] = windave_var_rs(vertvar_vary_rs, precip_occ)

% Minimize the vertical variance of RH averaged over the daily. Occurrence
% of daily precipitation (but not amount) is an optional input (used below)
% that is used as a guide in temporally smoothing estimates of r_surf,
% which otherwise fluctuate too much in response to synoptic weather
% variability

num_days = length(precip_occ);

% Estimate vertical variance relative to minimum (optional)
% min_val = repmat(min(vertvar_vary_rs),61,1);
% vertvar_vary_rs = vertvar_vary_rs - min_val;

% Set triangulate window size
wind_size = 21; % window size, wind_size (days) (has to be odd!)

% Initialize vector to store chosen r_surf index
RS = zeros(1,num_days)+NaN;

for DY = 1:(num_days) % DY counts trhough days
    
    % if there are any NaNs in vertvar_vary_rs(:,DY), that day was missing
    % forcing data (so keep RS == NaN for that day).
    if ~any(isnan(vertvar_vary_rs(:,DY)))
        
        % Calculate window size based on precipitation occurance:
        if precip_occ(DY) == 1 % If it rains, do not apply window
            wts=zeros(1,wind_size);
            wts(:,ceil(wind_size/2)) = ceil(wind_size/2);
        else % If it does not rain, apply 21-day triangular window
            wts=[1:floor(wind_size./2),ceil(wind_size./2),floor(wind_size./2):-1:1];
        end
        
        % Extract data (from vertvar_vary_rs) for 21-day period
        shift_wind = floor(length(wts)/2);
        arg_var = vertvar_vary_rs(:,min(max(1,DY-shift_wind:DY+shift_wind),num_days));
        
        % If there any NaNs in "arg_var", reweight wts
        wts=ones(61,1)*(wts./sum(wts));
        nanind=ones(size(wts)); % nanind indicates where there are NaNs
        nanind(isnan(arg_var))= NaN;
        
        % Weight variances
        ave_wt_wind = nansum(wts.*arg_var,2)./nansum(wts.*nanind,2);
        
        % Minimize the vertical varirance
        [~, index_RS] = min(ave_wt_wind);
        RS(DY) = index_RS; % save index of minimum
        
    end
    
end

end

