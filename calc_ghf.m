function ghf = calc_ghf(t,Isoil) 
% Estimated ground heat flux (GHF) from the time history of 
% estimated ground temperature using an analytical solution
% of the diffusion equation for heat transfer.
% For more details, see Appendix 1 of (Rigden and Salvucci, 2016, GCB)

disp('... runing ghf model')

B = (Isoil./sqrt(pi)).*2;

% Estimate ground temperature (tsg) from screen height air tempterature. 
% Note that C1 and C2 were calibrated using data from 20 AmeriFlux sites,
% thus, are not site specific.
C1 = 0.56; C2 = 2; 
t_daily = reshape(repmat(mean(t,1),[48,1,1]),48*365,1);
t_pert = repmat(t(:)-t_daily,2,1); % Here we add a leading year to minimize  
                                   % edge effects. If there are multiple 
                                   % years of data, use the previous year's  
                                   % temperature datat instead of 
                                   % repreating the same year twice.
temp_long_circ = circshift(t_pert,[C2,0]);
temp_long_fit = temp_long_circ.*C1;
tsg = temp_long_fit + repmat(t_daily,2,1);

% Smooth temperature with moving window (5-half hour span)
tsg_smth = smooth(tsg,5,'moving');

% Dummy variable
s = 1:1:(48*365*2);

% Initialize GHF
g_save = zeros(1,(length(s)-1));

for t=1:length(s)-1

    II = (t-(48*365*1):(t-1));
    i = II(II>0);

    DIF = (tsg_smth(i+1)-tsg_smth(i));
    SQRT = (1./(s(i+1)-s(i))').*(sqrt(t - s(i)) - sqrt(t - s(i+1)))'; 

    gg = DIF.*SQRT;
    g_save(t) = nansum(gg);

end

% Convert to correct units
% 1800 for units of "sqrt(dt)" and B = (Isoil./sqrt(pi)).*2;   
ghf_long = g_save*B./sqrt(1800);

% Delete leading year
ghf_long(1:(365*48-1)) = [];

% Reshape 
ghf = reshape(ghf_long,48,365);

end

