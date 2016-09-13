
num_yrs = 50;
s = 1:1:(48*365*(num_yrs+1));


% Initialize GHF
years_lag = 5;

G = zeros(1,(length(s)-1));

for kk=1:(48*365*years_lag+1)
    
    disp([kk (48*365*years_lag+1)])
    t=s(kk);
    
    gg = zeros(1,t);
 
    II = (kk-(48*365*years_lag):(kk-1));
    i = II(II>0);

    SQRT = (1./(s(i+1)-s(i))').*(sqrt(t - s(i)) - sqrt(t - s(i+1)))';     

end

save('SQRT_GHF2','SQRT')


