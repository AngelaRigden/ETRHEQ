function [balance_error, le, sh, rlu, L, chih, u_str, RIB] = ... 
    calc_energy_bal_error(ts, data, z_m, z_oh, z_ov, d, z_o, ...
    r_surf, den, u_pred_zm)
% function [energy_terms] = temp_bisectionG(ts, ts_1, data, z_m, z_m_wind, u_str, aaa, z_oh, z_ov, r_surf, den, iter, B,ghf)
% Given a surface temperature, time of day, and surface resistance, it will calculate the energy balance error.
global cp k g lv emis eps sb RIMAX

% S1 and S2 are used to reshape vectors and matrices
S1 = length(data(:,1)); % 48 if half hourly data
S2 = size(u_pred_zm,2); % 61 x num_frac

% Extract each variable
% Recall, data = [1-t, 2-q, 3-u, 4-p, 5-rsnet, 6-rld, 7-ghf];
t = data(:,1)*ones(1,S2);
q = data(:,2);
p = data(:,4);
rsnet = data(:,5);
rld   = data(:,6);
ghf   = data(:,7);

% To estimate (z-d)/L, use method oulined by:
% Blumel, K. (2000), An approximate analytical solution of flux-profile
% relationships for the atmospheric surface layer with different momentum
% and heat roughness lengths, Boundary-Layer Meteorol.

ym = ((z_m-d)./z_o)*ones(48,S2);
yh = (z_m-d)./z_oh; 
log_ym = log(ym);
log_yh = log(yh); 

% If bisection choose the actual air temperature, a numerical instability
% arises caused by introduction of zeros
temp_diff = t-ts;
temp_diff(abs(temp_diff) < (10^-10)) = (10^-10); % just make small & ~= 0

% Calculate flux richardson number
RIB = (g.*(z_m-d-z_o)./t).*temp_diff./(u_pred_zm.^2); % dry version
RIB = min(RIB,RIMAX);

ada=log_yh./log_ym;

% Stable conditions
% Eq. 27 in Blumel (2000)
y_pos1 = (ym./(ym-1)).*log_ym.*1.*RIB.*( ...
    (1-5.*1.*RIB) + ...
    ((ada/2)-1) + ...
    sqrt( ...
    1.*RIB.*(5.*(((yh-1)./(ym-1)).*(ym./yh))-5.*ada)+ ...
    0.25.*ada.^2) ...
    ).^-1;
y=y_pos1; 

% Unstable conditions
alpha=((yh-1)./(ym-1)).*(ym./yh);
CHSI=log_ym./log_yh;
betastar=(0.01857.*(log_ym-2.3026)+0.025).*((16./16).^10);
y_neg=(ym./(ym-1)).*1.*RIB.*CHSI.*log_ym.*sqrt(...
    (1-16.*RIB.*(CHSI.^2).*1.*alpha)./ ...
    (1-16.*RIB.*(CHSI).*1)).* ...
    (1./(1-betastar.*RIB.*1));
SOL4=find(RIB<0);
y(SOL4)=y_neg(SOL4);

% now, with (z-d)/L calculated, find the psih and psim functions
L = (z_m-d)./y;

xi_vh = ((z_m-d)./L);   % dummy-variable for vapor and heat
xi_v  = ((z_ov)./L);     % dummy-variable for vapor
xi_h  = ((z_oh)./L);     % dummy-variable for heat
xi_m1 = ((z_m-d)./L);   % dummy-variable for momentum
xi_m2 = (z_o./L);       % dummy-variable for momentum

% Estimate stability corrections
stbl_vh = zeros(S1,S2)+NaN;
stbl_v  = zeros(S1,S2)+NaN;
stbl_h  = zeros(S1,S2)+NaN;
stbl_m1 = zeros(S1,S2)+NaN;
stbl_m2 = zeros(S1,S2)+NaN;

phi_vh = zeros(S1,S2)+NaN;
phi_v  = zeros(S1,S2)+NaN;
phi_h  = zeros(S1,S2)+NaN;
phi_m1 = zeros(S1,S2)+NaN;
phi_m2 = zeros(S1,S2)+NaN;

% VH
II=find(xi_vh<=0);
phi_vh(II)=(1-16.*xi_vh(II)).^0.5;
stbl_vh(II) = 2.*log((1+phi_vh(II))./2);
stbl_vh(xi_vh>0) = -5.*xi_vh(xi_vh>0);

% V
II=find(xi_v<=0);
phi_v(II)=(1-16.*xi_v(II)).^0.5;
stbl_v(II) = 2.*log((1+phi_v(II))./2);
stbl_v(xi_v>0) = -5.*xi_v(xi_v>0);

% H
II=find(xi_h<=0);
phi_h(II)=(1-16.*xi_h(II)).^0.5;
stbl_h(II) = 2.*log((1+phi_h(II))./2);
stbl_h(xi_h>0) = -5.*xi_h(xi_h>0);

% M1
II=find(xi_m1<=0);
phi_m1(II)=(1-16.*xi_m1(II)).^0.25;
stbl_m1(II) = 2.*log((1+phi_m1(II))./2) + log((1+phi_m1(II).^2)./2) - 2.*atan(phi_m1(II)) + pi./2;
stbl_m1(xi_m1>0) = -5.*xi_m1(xi_m1>0);

% M2
II=find(xi_m2<=0);
phi_m2(II)=(1-16.*xi_m2(II)).^0.25;
stbl_m2(II) = 2.*log((1+phi_m2(II))./2) + log((1+phi_m2(II).^2)./2) - 2.*atan(phi_m2(II)) + pi./2;
stbl_m2(xi_m2>0) = -5.*xi_m2(xi_m2>0);

chiv = (log_yh - stbl_vh + stbl_v);
chih = (log_yh - stbl_vh + stbl_h); 
chim = (log_ym - stbl_m1 + stbl_m2);

% Near the surface, the potential temperature is ~= to the temperature.
ts_pot = ts;

% Estimate friction velocity, u_str (m/s)
u_str = u_pred_zm.*k./(chim);

% Calculate sensible heat
sh = ((u_str.*(den*ones(1,S2)).*cp.*k).*(ts_pot-t))./chih;

% Calculate saturation vapor pressure using the integrated
% clausius clayperon relation
es_sat = 611.2.*exp((17.67.*(ts_pot-273.15))./(ts_pot-29.65)); 

% To keep from boiling, *rarely* kicks in (but since we are running for a
% large range of r_surf values, sometimes odd solutions occur)
es_sat = min(es_sat,0.95*(p*ones(1,S2)));

% Calculate saturated specific humidity at the surface
qs_sat = (eps.*es_sat)./((p*ones(1,S2))-((1-eps).*es_sat)); 

% Calculate et
et_final = (qs_sat-(q*ones(1,S2)))./((r_surf./(den*ones(1,S2))) + (chiv./(k.*u_str.*(den*ones(1,S2)))) );

% Calculate latent heat
le = lv.*et_final;

% Calculate sensible heat
rlu = emis.*sb.*(ts_pot.*ts_pot.*ts_pot.*ts_pot);

% Calculate energy balance error
rlnet = bsxfun(@minus,rld,rlu);
rnet = bsxfun(@plus,rsnet,rlnet);
balance_error = rnet - le - ghf*ones(1,S2) - sh;

end
