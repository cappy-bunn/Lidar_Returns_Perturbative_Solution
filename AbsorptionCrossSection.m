% Absorption Cross Section Calculation
function [a,b] = AbsorptionCrossSection(wl_laser,wl_absorption,T,P,r)
% Where wl_absorption is the online wavelength in nm
% wl_laser is the incident light wavelength (nm)
% T is the temperature profile
% P is the pressure profile
% r is the heights vector
% This returns a, which is the cross section of absorbers at the online wavelength's frequency
% and b, which is a matrix that gives the cross section of absorbers at any frequency in the range.
    
c = 3.00e8;          % Speed of light [m/s]
kb = 1.38065e-23;    % Boltzmann constant [J/K]
h = 6.626e-34;       % Planck constant [J s]
m = 5.314e-26;       % Mass of O2 molecule [kg]

wn0 = 1e9/wl_laser;                 % Wavenumber [m^-1]
fwidth = 20;                        % Range to scan [GHz]
fstepnum = 301;                     % Number of steps
fstart = c*wn0-fwidth*1e9/2;        % Beginning of freq range [s^-1]
fstep = fwidth*1e9/(fstepnum-1);    % Step width [s^-1]
fstop = c*wn0+fwidth*1e9/2;         % End of freq range [s^-1]

nu0 = 1/(wl_absorption*10^(-9)*100);% Wavenumber [cm^-1] corresponding to center of absorption line
S0 = 1.08e-25;                      % (HITRAN) line strength [cm/molecule]
t0 = 293.15;                        % Reference temperature [K]
p0 = 1;                             % Reference pressure [atm]
ccm = 3e10;                         % Speed of light [cm/s]
E = 1248.204;                       % Energy above the g.s. of the absorption line's lower level [cm^-1]
gl0 = 0.0362;                       % Halfwidth [cm^-1] (gl0*c = 1.085 GHz)


j=1;
for freq = fstart:fstep:fstop
    f(j) = (freq-wn0*c)/1e9;        % Distance from laser wavenumber [GHz]
    lambda(j) = c/freq;             % Wavelengths within range [m]
    j=j+1;
end

k = length(r);

for j=1:1:k
    t = T(j);
    pa = P(j);
    a = S0*(t0/t);
    b = (1-exp(-h*ccm*nu0/(kb*t)))/(1-exp(-h*ccm*nu0/(kb*t0)));
    d = exp(h*ccm*E*(1/t0-1/t)/kb);
    st = a*b*d;                                         % Nehrir Eq. 2.15
    
    gammaL = gl0*(pa/p0)*(t0/t)^0.71;                   % Nehrir Eq. 2.16
    gammaD = (nu0/ccm)*(2*kb*t*log(2)/m)^0.5*100;       % Nehrir Eq. 2.21
    i=1;
    
    for freq = fstart:fstep:fstop
        nu = freq/c/100;
        
        % Calculate the convolution integral
        x = (nu-nu0)*0.8325546/gammaD;
        y = gammaL*0.8325546/gammaD;
        ttt = -3:0.0001:3;
        ft = exp(-ttt.^2)./(y.*y+(x-ttt).^2);
        convint = trapz(ttt,ft);
        sigmaV(j,i) = st*0.12448*gammaL/(gammaD^2)*convint;
        
        i=i+1;
    end
    
end

a = (sigmaV(:,151))';

b = sigmaV;
end