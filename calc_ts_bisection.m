function [TS, LE, SH, RLU, L, CHIH, U_STR, RIB] = ...
    calc_ts_bisection(data,...
    r_surf, den, ...
    z_m, z_o, z_ov, z_oh, d, ...
    u_pred_zm)

% The function "calc_ts_bisection" uses the method of bisection to find the
% roots of the energy blance equations (where the equation balances to
% zero).

% The bisection algorithm runs until:
% (1) energy balance error for each half hour is < enbal_error_threshold
% (2) maximum number of iterations is reached (when there are errors in the
%     input measurements, this may occur). An error will be thrown and the
%     data should be looked at.
enbal_error_threshold = 1;   % energy balance error threshold (W/m^2)
max_number_iterations = 200; % maximum number of iterations in bisection

% Initial bounding surface temperatures used in bisection
ts1 =  100*ones(size(r_surf)); % lower ts bound
ts2 =  400*ones(size(r_surf)); % upper ts bound

% Calculate energy balance error
bal_err1 = calc_energy_bal_error(ts1, data, z_m, z_oh, z_ov, d, z_o, ...
    r_surf, den, u_pred_zm); % ts1 & ts2 need to be NN x RSN arrays.
bal_err2 = calc_energy_bal_error(ts2, data, z_m, z_oh, z_ov, d, z_o, ...
    r_surf, den, u_pred_zm);

% Determine sign of energy balance error
sign_balerr1 = sign(bal_err1);
sign_balerr2 = sign(bal_err2);

% Final temperature array
TS = zeros(size(r_surf));

iter = 0;
Z = 1;
while (iter < max_number_iterations && Z ~= 0)
    
    % Find the midpoint of the bounding temperatures
    ts_mid = (ts1 + ts2)/2;
    bal_errmid = calc_energy_bal_error(ts_mid, data, z_m, z_oh, z_ov, d, z_o, ...
        r_surf, den, u_pred_zm);
    
    % Check and see if mid-point is within error.
    II = find(abs(bal_errmid) < enbal_error_threshold);
    TS(II) = ts_mid(II);
    
    II = find(sign_balerr2 == sign(bal_errmid));
    ts2(II)  = ts_mid(II);
    bal_err2(II) = bal_errmid(II);
    
    II = find(sign_balerr1 == sign(bal_errmid));
    ts1(II) = ts_mid(II);
    bal_err1(II) = bal_errmid(II);
    
    Z = length(find(TS(:)==0));
    iter=iter+1;
end

if (iter == max_number_iterations)
    error('Surface temperature bisection did not converge')
end

ts_mid = (ts1 + ts2)/2;

[BAL_ERR, LE, SH, RLU, L, CHIH, U_STR, RIB] = ...
    calc_energy_bal_error(ts_mid, data, z_m, z_oh, z_ov, d, z_o, ...
    r_surf, den, u_pred_zm);
TS = ts_mid;

end
