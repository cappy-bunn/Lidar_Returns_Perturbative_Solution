clear all; %#ok<CLALL>
close all;
clc;

%-------- Definitions -------

global wl_on wl_off Ts Ps bg bl blw sa sm c n0 kb h m Epulse tpulse...
    rangebin r rkm n fwidth wn0_on f0_on f1_on f2_on f_on f wn1_on...
    wn2_on wn_on_step wn_on wn0_off f0_off f1_off f2_off wn1_off wn2_off...
    wn_off_step wn_off;

% Enter online and offline wavelengths [nm]:
wl_on = 769.2332;
wl_off = 769.1332;

Ts = 293.15;    % Surface temperature [K]
Ps = 0.84;      % Surface pressure [atm] at 4793 ft (Bozeman)
gamma = -7;     % Adiabatic lapse rate [K/km]

bg = 0.3e-6;    % Backscatter at ground
bl = 3;         % Boundary layer height [km]
blw = .75;      % Boundary layer falloff width [km]
sa = 50;        % Aerosol lidar ratio
sm = 8.3776;    % Molecular lidar ratio (8/3)*pi

c = 3.00e8;          % Speed of light [m/s]
n0 = 2.687e25;       % Loschmidt constant [1/m^3] at 0 C and 1 atm
kb = 1.38065e-23;    % Boltzmann constant [J/K]
h = 6.626e-34;       % Planck constant [J s]
m = 5.314e-26;       % Mass of O2 molecule [kg]

Epulse = 1e-5;          % Pulse energy [J]
tpulse = 1e-6;          % Pulse duration [s]
rangebin = c*tpulse/2;  % Range bin size [m]
r = 300:rangebin:6000;  % Heights [m]
rkm = r./1000;          % Heights [km]

n = 301;                % Number of steps in range
fwidth = 20;            % Range to scan [GHz]

% ONLINE
wn0_on = 1e9/wl_on;             % Wavenumber of ONLINE laser [m^-1]
f0_on = wn0_on*c/(1e9);         % Frequency of ONLINE laser [GHz]
f1_on = f0_on-fwidth/2;         % Start of frequency range [GHz]
f2_on = f0_on+fwidth/2;         % End of frequency range [GHz]
f_on = f1_on:((f2_on-f1_on)/(n-1)):f2_on; % Frequency range [GHz]
f = f_on-f0_on;                 % Frequency range wrt online wl [GHz]
wn1_on = (f1_on*1e9)/c;         % Start of wavenumber range [m^-1]
wn2_on = (f2_on*1e9)/c;         % End of wavenumber range [m^-1]
wn_on_step = (wn2_on-wn1_on)/(n-1); % Step size
wn_on = wn1_on:wn_on_step:wn2_on; % Wavenumbers [m^-1]


% OFFLINE
wn0_off = 1e9/wl_off;           % Wavenumber of OFFLINE laser [m^-1]
f0_off = wn0_off*c/(1e9);       % Frequency of OFFLINE laser [GHz]
f1_off = f0_off-fwidth/2;       % Start of frequency range [GHz]
f2_off = f0_off+fwidth/2;       % End of frequency range [GHz]
wn1_off = (f1_off*1e9)/c;       % Start of wavenumber range [m^-1]
wn2_off = (f2_off*1e9)/c;       % End of wavenumber range [m^-1]
wn_off_step = (wn2_off-wn1_off)/(n-1);  % Step size
wn_off = wn1_off:wn_off_step:wn2_off; % Wavenumbers [m^-1]


%% --------- Profiles ---------
T = Ts+gamma.*rkm;                      % Temperature [K]
P = Ps*(Ts./T).^(-5.2199);              % Pressure [atm]
Pkpa = 101.325*P;                       % Pressure [kPa]
nO2 = 0.2096*n0*T.*P/(273.15*1);        % O2 number density [m^-3]

betam1 = 374280*Pkpa./(T*wl_on^4);      % ONLINE Molecular Backscatter
extm1 = sm*betam1;                      % ONLINE Molecular Extinction (Casey ~Eq. 45)

betam2 = 374280*Pkpa./(T*wl_off^4);     % OFFLINE Molecular Backscatter
extm2 = sm*betam2;                      % OFFLINE Molecular Extinction

% Aerosol Backscatter
k = size(r);
k = k(2);
for j=1:1:k
    range = rkm(j);
    betaa(j) = bg/100;
    if range<bl+blw/2
        betaa(j) = bg*(1-(bl-blw/2)/10)*(1-(range-(bl-blw/2))/blw);
    end
    if range<bl-blw/2
        betaa(j) = bg*(1-range/10);
    end
end
exta = sa*betaa;                      % Aerosol Extinction (Casey ~Eq. 45)

% ONLINE Doppler Broadened Backscatter Lineshape
D1_unnorm = DopplerBroadBeta(wn0_on,wn_on,m,T);
D1_norm_func = trapz(D1_unnorm.*wn_on_step);
D1 = D1_unnorm./D1_norm_func./2;                   % Normalized lineshape

% OFFLINE Doppler Broadened Backscatter Lineshape
D2_unnorm = DopplerBroadBeta(wn0_off,wn_off,m,T);
D2_norm_func = trapz(D2_unnorm.*wn_off_step);
D2 = D2_unnorm./D2_norm_func./2;                   % Normalized lineshape
%% -------- Absorption Cross Section Calculation ---------
[cs1,sigmaV1] = AbsorptionCrossSection(wl_on,wl_on,T,P,r);
[cs2,sigmaV2] = AbsorptionCrossSection(wl_off,wl_on,T,P,r);
alpha1 = cs1.*nO2/(100^2);                      % absorption coefficient for online wavelength
alpha2 = cs2.*nO2/(100^2);                      % absorption coefficient for offline wavelength
alpham1 = sigmaV1.*repmat(nO2',1,301)/(100^2);  % repeat array nO2 to create matrix same size as sigmaV1
alpham2 = sigmaV2.*repmat(nO2',1,301)/(100^2);  % repeat array nO2 to create matrix same size as sigmaV2
% alpham is absorption coefficient for given height and frequency (within our range)

%% -------- Calculate Returns ---------
% First find optical depth from absorption coefficient (alpham).
optdepth1 = OpticalDepth(alpham1,rangebin);  % Optical depth for online
optdepth2 = OpticalDepth(alpham2,rangebin);  % Optical depth for offline

% Next find molecular absorption: Tm = exp[-int(alpha)dr from r0 to r]
Tm1 = exp(-optdepth1);       % Fraction of online light reaching a range r
Tm2 = exp(-optdepth2);       % Fraction of offline light reaching a range r

% Find atmospheric transmission: Ta = exp[-int(total extinction)dr from r0 to r
Ta1 = AtmosTransmission(extm1,exta);
Ta2 = AtmosTransmission(extm2,exta);

% This was for a non-delta laser lineshape.
% FWHM of diode laser is about 1.5 nm. Translates to 0.0076 GHz. 
% FWHM = 2.355*stddev
%   stddev = 0.00322718;            % standard deviation [GHz]
%   hi = exp((-f.^2)/(2*stddev^2)); % Laser spectral line shape

%   int_Tu = hi'.*Tm';              % Integrad of Tu integral. Transposed so 'trapz' integrates over frequency
%   Tu = trapz(int_Tu);             % Transmittance going up thru atmosphere


% For a delta laser lineshape hi = delta(nu-nu0),
% Tu is int(hi Tm) over frequency, but delta function picks out nu0:
Tu1 = Tm1(:,151)';                % Fraction of online light reaching range r w/ lineshape
Tu2 = Tm2(:,151)';                % Fraction of offline light reaching range r w/ lineshape

E_on = f0_on*(1e9)*h;           % Energy of online photon
N01 = Epulse/E_on;              % Number of photons initially emitted online if pulse is 10 uJ
E_off = f0_off*(1e9)*h;         % Energy of offline photon
N02 = Epulse/E_off;             % Number of photons initially emitted offline if pulse is 10 uJ

Nu1 = N01*Ta1.*Tu1;              % Number of online photons that reach the scatterer
Nu2 = N02*Ta2.*Tu2;              % Number of offline photons that reach the scatterer

betat1 = betaa+betam1;          % Total online backscatter
betat2 = betaa+betam2;          % Total offline backscatter

aero_back1 = repmat((betaa./betat1),301,1);     % Aerosol reflects delta function lineshape
aero_back1(1:150,:)=0;
aero_back1(152:301,:)=0;

gi1 = aero_back1+(betam1./betat1).*D1;  % Total return lineshape online

aero_back2 = repmat((betaa./betat2),301,1);     % Aerosol reflects delta function lineshape
aero_back2(1:150,:)=0;
aero_back2(152:301,:)=0;

gi2 = aero_back2+(betam2./betat2).*D2;  % Total return lineshape offline

Ns1 = Nu1.*betat1.*gi1;             % Number of online photons that backscatter
Ns2 = Nu2.*betat2.*gi2;             % Number of offline photons that backscatter

A = 0.25;                       % Area of telescope entrance pupil [m^2]
Nt1 = Nu1.*betat1.*gi1.*Ta1.*Tm1'.*(A./r.^2)*rangebin;    %Number of online photons seen by the telescope
Nt2 = Nu2.*betat2.*gi2.*Ta2.*Tm2'.*(A./r.^2)*rangebin;    %Number of offline photons seen by the telescope

overlap = 1;                    % Overlap function
epsilon = 0.25;                 % Receiver transmission*detector efficiency*max etalon transmission

% Etalon Transmission calculations
wl_start = wl_on-0.01972;           % Beginning of wavelength range [nm] (equivalent to 10 GHz away from online)
wl_end = wl_on+0.01972;             % End of wavelength range [nm]
wl_step = (wl_end-wl_start)/(n-1);  % Wavelength step size (made so the range has 301 points)
lambda = wl_start:wl_step:wl_end;   % Wavelength range

R = 0.92956;                            % Etalon reflectivity
FSR = wl_off-wl_on;                     % Free Spectral Range (FSR=c/(2*nL)) difference between online and offline wl [nm]
theta = (2*pi*(lambda-wl_on))/FSR;                      % Round trip phase accumulation
etalonT_unnorm = (1-R)^2./(pi*(1+R^2-2*R*cos(theta)));  % Etalon transmission (before normalization)
etalon_norm_factor = etalonT_unnorm(:,151);             % Normalization factor
etalonT = etalonT_unnorm/etalon_norm_factor;            % Normalized etalon transmission
%etalonT = 1;
Nd1 = Nt1'.*overlap*epsilon.*etalonT;  % Number of photons seen by the detector
Nd2 = Nt2'.*overlap*epsilon.*etalonT;  % Number of photons seen by the detector

int_Td1 = gi1.*Tm1'.*etalonT';
Td1 = trapz(int_Td1);             % Transmittance down thru atmosphere

int_Td2 = gi2.*Tm2'.*etalonT';
Td2 = trapz(int_Td2);             % Transmittance down thru atmosphere

Nr1 = Nu1.*betat1.*Ta1.*overlap*epsilon.*(A./r.^2)*rangebin.*Td1;  % Number of photons seen by the receiver
Nr2 = Nu2.*betat2.*Ta2.*overlap*epsilon.*(A./r.^2)*rangebin.*Td2;  % Number of photons seen by the receiver

%semilogx(Nr1,rkm,Nr2,rkm)               % Log plot of Nr (# of photons) vs height (km)

%% The DIAL Equation
% Find the derivative wrt r of gi1:
num_col = size(gi1,2);
ii=1:num_col-2;
jj=3:num_col;
for ix = 1:length(ii)
    ii(ix);
    jj(ix);
    dgi1_dr_new(:,ii) = (gi1(:,jj)-gi1(:,ii))/(2*rangebin);
end

% Repeat the first and last columns to maintain original matrix dimensions:
dgi1_dr = [dgi1_dr_new(:,1) dgi1_dr_new dgi1_dr_new(:,size(dgi1_dr_new,2))];

int_numG1 = trapz(dgi1_dr.*etalonT'.*Tm1');
int_denG1 = trapz(gi1.*etalonT'.*Tm1');
G1 = int_numG1./int_denG1;         % Correction factor

% Find the derivative wrt r of gi2:
num_col2 = size(gi2,2);
ii=1:num_col2-2;
jj=3:num_col2;
for ix = 1:length(ii)
    ii(ix);
    jj(ix);
    dgi2_dr_new(:,ii) = (gi2(:,jj)-gi2(:,ii))/(2*rangebin);
end

% Repeat the first and last columns to maintain original matrix dimensions:
dgi2_dr = [dgi2_dr_new(:,1) dgi2_dr_new dgi2_dr_new(:,size(dgi2_dr_new,2))];

int_numG2 = trapz(dgi2_dr.*etalonT'.*Tm2');
int_denG2 = trapz(lambda'.*gi2.*etalonT'.*Tm2');
G2 = int_numG2./int_denG2;         % Correction factor

int_num_alphaD1 = trapz(gi1.*etalonT'.*alpham1'.*Tm1');
int_den_alphaD1 = trapz(gi1.*etalonT'.*Tm1');
alphaD1 = int_num_alphaD1./int_den_alphaD1;

int_num_alphaD2 = trapz(gi2.*etalonT'.*alpham2'.*Tm2');
int_den_alphaD2 = trapz(gi2.*etalonT'.*Tm2');
alphaD2 = int_num_alphaD2./int_den_alphaD2;

int_num_alphaU1 = alpham1(:,151).*Tm1(:,151);
int_den_alphaU1 = Tm1(:,151);
alphaU1 = int_num_alphaU1./int_den_alphaU1;

int_num_alphaU2 = alpham2(:,151).*Tm2(:,151);
int_den_alphaU2 = Tm2(:,151);
alphaU2 = int_num_alphaU2./int_den_alphaU2;

DIAL_exact = -alphaU1'+alphaU2'-alphaD1+alphaD2+G1-G2;

% Where alphaD2 approx = alphaU2, so create new term alphag2 = alphaD2 = alphaU2
alphag2 = alphaD2;

%% -------- Model G1 ---------
% G1 can be written as the sum of a modeled G1 and a small perturbation.
% First find modeled G1 (G1m).
gamma_model = -6;
G1m = ModeledG1(gamma_model);

%% -------- Find W(nu,r) ---------
alphag1 = alpham1(:,151);                           % New DIAL Notes sec. 3
f_ls = alpham1./alphag1;                            % Normalized line shape
                                                    % New DIAL Notes sec. 4

int_num_W = trapz(gi1.*etalonT'.*(f_ls'-1).*Tm1'*wl_step);
int_den_W = trapz(gi1.*etalonT'.*Tm1'*wl_step);

W = int_num_W./int_den_W;                           % New DIAL Notes sec. 4

%% --------- Perturbative Solution -----------
% ZEROTH
r2 = 2:size(alphag2,2);
r1 = r2-1;
ln_returns = (log((Nr1(r2).*Nr2(r1))./(Nr1(r1).*Nr2(r2))))./rangebin;
alpha0_b = alphag2(:,1:size(ln_returns,2))-0.5*ln_returns;
alpha0 = [alpha0_b alpha0_b(:,size(alpha0_b,2))];                 % Repeat last row to maintain original matrix size
% plot(alpha0,rkm,alphag1,rkm)


% FIRST
% Calculate G1 and W using alpha0 where Tm = exp(-int(alpha0)dr from r0 to r)
optdepth_approx = OpticalDepth(alpha0'.*f_ls,rangebin);
Tm_approx = exp(-optdepth_approx);

% Recalculate G1 using Tm_approx
int_numG1_approx = trapz(dgi1_dr.*etalonT'.*Tm_approx');
int_denG1_approx = trapz(gi1.*etalonT'.*Tm_approx');
G1_approx = int_numG1_approx./int_denG1_approx;

% Recalculate W using Tm_approx
int_numW_approx = trapz(gi1.*etalonT'.*(f_ls'-1).*Tm_approx'*wl_step);
int_denW_approx = trapz(gi1.*etalonT'.*Tm_approx'*wl_step);
W_approx = int_numW_approx./int_denW_approx;  

% First order correction
delta_alpha1 = 0.5*(G1_approx-alpha0.*W_approx);

rkms = 0.375:0.15:6.1;
% plot(alpha0,rkms,alphag1',rkm)
% plot(alpha0+delta_alpha1,rkms,alphag1',rkm)

rkm2 = rkms(:,1:38);
% plot(ln_returns,rkm2,DIAL_exact,rkm)

% SECOND
% Calculate next order term of Tm (Tm_approx_1) from delta_alpha1, then find delta_W_approx and
% delta_G1_approx.
optdepth_approx_1 = OpticalDepth(delta_alpha1'.*f_ls,rangebin);
Tm_approx_1 = optdepth_approx_1.*Tm_approx;

int_numW_approx_1 = trapz(gi1.*etalonT'.*(f_ls'-1).*Tm_approx_1'*wl_step);
delta_W_approx = -int_numW_approx_1./int_denW_approx;

int_numG1_approx_1 = trapz(dgi1_dr.*etalonT'.*Tm_approx_1');
delta_G1_approx = -int_numG1_approx_1./int_denG1_approx;

delta_alpha2 = 0.5*(delta_G1_approx - delta_alpha1.*W_approx - alpha0.*delta_W_approx);

% plot(alpha0+delta_alpha1+delta_alpha2,rkms,alphag1',rkm)
% plot(alpha0+delta_alpha1,rkms,alphag1',rkm)
