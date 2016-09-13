function r_ld = calc_rld(p, t, q, z_elev)
% This function estimates downwelling longwave radiation as a function of
% screen height air temperature, vapor pressure, and elevation 

global eps

% Calculate vapor pressure (e, Pa) at surface from specific humidity, q,
% and pressure, p:
e = (q.*p)./(q.*(1-eps)+eps);

% Abramowitz method: Abramowitz, G., L. Pouyanné, and H. Ajami (2012), On
% the information content of surface meteorology for downward atmospheric
% long?wave radiation synthesis, Geophys. Res. Lett., 39(4),
% doi:10.1029/2011GL050726. 
% LW = 0.031.*e + 2.84.*t - 522.5;

% Revised Abramowitz method, described in Appendix 1 of (Rigden and
% Salvucci, 2016, GCB)
r_ld = 0.045.*e + 2.18.*t - 0.0069*z_elev - 351.2;

end