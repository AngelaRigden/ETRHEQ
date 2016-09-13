function PSIsave = calc_stability_corr(z_val,d,z_o,L)

% Calculate stability correction
xi_m1 = ((z_val-d)./L); % dummy-variable for momentum
xi_m2 = (z_o./L);       % dummy-variable for momentum

stbl_m1 = zeros(size(L))+NaN;
stbl_m2 = zeros(size(L))+NaN;
phi_m1  = zeros(size(L))+NaN;
phi_m2  = zeros(size(L))+NaN;

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

PSIsave = stbl_m1-stbl_m2;

end

