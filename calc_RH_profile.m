function RHz = calc_RH_profile(data, u_str_new, ts, h, le, L, z, ...
    r_surf, t_ref, den, z_oh, z_ov, d)
%This function calculated the relative humidity profile.
% OUTPUTS = RH, theta_z, tz, qz, qz_sat

global r_d cp k g eps lv 

% size(z) = 48 x 61 x length(z_levels)
% number of vertical (z) levels, num_levs
num_levs = size(z, 3);

r_surf = r_surf(:,:,ones(1,1,num_levs));

% Load surface pressure
p = data(:,4)*ones(1,61);

u_str = u_str_new(:,:,ones(1,1,num_levs));
p = p(:,:,ones(1,1,num_levs));

ts = ts(:,:,ones(1,1,num_levs));
le = le(:,:,ones(1,1,num_levs));
h = h(:,:,ones(1,1,num_levs));
L = L(:,:,ones(1,1,num_levs));

z_oh = z_oh(:,:,ones(1,1,num_levs));
z_ov = z_ov(:,:,ones(1,1,num_levs));
den = den(:,:,ones(1,1,num_levs));

% Calculate saturateion vapor pressure using the integrated
% clausius clayperon relation at the surface
es_sat = 611.2.*exp((17.67.*(ts-273.15))./(ts-29.65)); %S10

% Calculate saturated specific humidity at the surface
qs_sat = (eps.*es_sat)./(p-((1-eps).*es_sat)); %S9

% Calculate qs
qs = qs_sat - ((r_surf.*(le./lv))./den);

% Account for condensation
con = find(le<0);
qs_end = qs(:,61,:);
qs_end = qs_end(:,:,ones(1,61,num_levs));
qs(con) = qs_end(con);


% Potential temperature at the surface == ts
theta_s = ts; % because z = 0.

% Determine stability of atmosphere @ each z. 

% Calculate three dimensionless hieghts
xi_vh = (z-d)./L;  % dummy-variable for vapor and temp
xi_v = (z_ov)./L;    % dummy-variable for vapor
xi_h = (z_oh)./L;    % dummy-variable for temp

% Determine stability based on dimensionless hieghts (above) and equation(s) S4
    
% Stabilitly...
stbl_vh=zeros(48,61,num_levs)+NaN;
stbl_v=zeros(48,61,num_levs)+NaN;
stbl_h=zeros(48,61,num_levs)+NaN;
phi_vh=zeros(48,61,num_levs)+NaN;
phi_v=zeros(48,61,num_levs)+NaN;
phi_h=zeros(48,61,num_levs)+NaN;

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

chiv = (log((z-d)./z_ov) - stbl_vh + stbl_v);
chih = (log((z-d)./z_oh) - stbl_vh + stbl_h);

% Calculate q at height z
qz = qs - (((le./lv)./(k.*u_str.*den)).*chiv);

% Calculate potential temperature at each height z
theta_z = theta_s - ((h./(k.*u_str.*den.*cp)).*chih);

% Calculate T from potential temperatures at each height z
tz = theta_z.*((exp((-g.*z)./(r_d.*t_ref))).^(r_d./cp));

% Calculate the saturation vapor pressure at each height z
ez_sat = 611.2.*exp((17.67.*(tz-273.15))./(tz-29.65));

% Calculate the pressure at height z
p_z = p.*exp((-g.*z)./(r_d.*t_ref));

% Calculate the specific humidity at each height z
qz_sat = (eps.*ez_sat)./(p_z-((1-eps).*ez_sat));

% calculate relative humidity at each height
RHz = qz./qz_sat;

end

